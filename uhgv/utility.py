##
## misc functions
##

import csv
import gzip
import logging
import multiprocessing as mp
import os
import platform
import resource
import shutil
import signal
import subprocess as sp
import sys
import time
from collections import defaultdict
from operator import itemgetter

import Bio.SeqIO
import psutil


def get_logger(quiet):
    logger = logging.getLogger(__name__)
    if not quiet:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    formatter = logging.Formatter(fmt="%(message)s")
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(stream_handler)
    return logger


def max_mem_usage():
    """Return max mem usage (GB) of self and child processes"""
    max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if platform.system() == "Linux":
        return (max_mem_self + max_mem_child) / float(1e6)
    else:
        return (max_mem_self + max_mem_child) / float(1e9)


##
## misc utilities
##


def mean(values):
    return sum(values) / len(values)


def split_dmnd(inpath, outdir, num_splits, ext=""):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    total_size = os.stat(inpath).st_size
    split_size = int(total_size / num_splits)

    last_id = None
    split_num = 1
    cursize = 0
    outfile = open(os.path.join(outdir, str(split_num)), "w")
    infile = open(inpath)
    for l in infile:
        outfile.write(l)
        cursize += len(l)
        last_id = l.split()[0].rsplit("_", 1)[0]
    for l in infile:
        cur_id = l.split()[0].rsplit("_", 1)[0]
        if cursize > split_size and cur_id != last_id:
            split_num += 1
            cursize = 0
            outfile = open(os.path.join(outdir, str(split_num)), "w")
        outfile.write(l)
        cursize += len(l)
        last_id = cur_id
    outfile.close()   
    infile.close()


def split_fasta(inpath, outdir, num_splits, ext):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    total_size = 0
    handle = gzip.open(inpath, "rt") if inpath.endswith(".gz") else open(inpath)
    for r in Bio.SeqIO.parse(handle, "fasta"):
        total_size += len(r.seq)

    split_size = int(total_size / num_splits)

    split_num = 1
    cursize = 0
    out = open(os.path.join(outdir, str(split_num)) + ext, "w")

    handle = gzip.open(inpath, "rt") if inpath.endswith(".gz") else open(inpath)
    for r in Bio.SeqIO.parse(handle, "fasta"):
        if cursize > split_size:
            split_num += 1
            cursize = 0
            out = open(os.path.join(outdir, str(split_num)) + ext, "w")
        out.write(">" + r.id + "\n" + str(r.seq) + "\n")
        cursize += len(r.seq)

    out.close()


def parallel_prodigal(tmpdir, input, output, threads, cleanup):

    tmpdir_prodigal = os.path.join(tmpdir, "prodigal")

    # split input DNA seqs into chunks
    split_fasta(inpath=input, outdir=tmpdir_prodigal, num_splits=threads, ext=".fna")

    # list shell commands
    commands = []
    for file in os.listdir(tmpdir_prodigal):
        if not file.endswith(".fna"):
            continue
        tmpin = os.path.join(tmpdir_prodigal, file)
        tmpout = os.path.join(tmpdir_prodigal, file.rsplit(".", 1)[0] + ".faa")
        cmd = "prodigal-gv -p meta "
        cmd += f"-i {tmpin} "
        cmd += f"-a {tmpout} "
        cmd += f"1> /dev/null "
        cmd += f"2> {tmpout}.log"
        commands.append([cmd])

    # run shell commands in parallel
    return_codes = parallel(run_shell, commands, threads)
    if sum(return_codes) != 0:
        msg = "\nError: One or more prodigal-gv tasks failed to run\n"
        msg += f"See logs for details: {tmpdir_prodigal}/*.log"
        sys.exit(msg)

    # combine outputs
    with open(output, "w") as out:
        for file in os.listdir(tmpdir_prodigal):
            if not file.endswith(".faa"):
                continue
            for line in open(os.path.join(tmpdir_prodigal, file)):
                out.write(line)

    # clean up
    if cleanup:
        shutil.rmtree(tmpdir_prodigal)


##
## code for parallelization
##

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def terminate_tree(pid, including_parent=True):
    parent = psutil.Process(pid)
    for child in parent.children(recursive=True):
        child.terminate()
    if including_parent:
        parent.terminate()


def run_shell(cmd):
    p = sp.Popen(cmd, shell=True)
    return p.wait()

        
def parallel(function, arguments_list, threads):
    pool = mp.Pool(threads, init_worker)
    try:
        results = []
        for arguments in arguments_list:
            p = pool.apply_async(function, args=arguments)
            results.append(p)
        pool.close()
        while True:
            if all(r.ready() for r in results):
                return [r.get() for r in results]
            time.sleep(1)
    except KeyboardInterrupt:
        pid = os.getpid()
        terminate_tree(pid)


##
## code to calculate ANI
##


def parse_blast(handle):
    for line in handle:
        r = line.split()
        yield {
            "qname": r[0],
            "tname": r[1],
            "pid": float(r[2]),
            "len": float(r[3]),
            "qcoords": sorted([int(r[6]), int(r[7])]),
            "tcoords": sorted([int(r[8]), int(r[9])]),
            "qlen": float(r[-2]),
            "tlen": float(r[-1]),
            "evalue": float(r[-4]),
        }


def yield_alignment_blocks(handle):
    # init block with 1st record
    key, alns = None, None
    for aln in parse_blast(handle):
        key = (aln["qname"], aln["tname"])
        alns = [aln]
        break
    # loop over remaining records
    for aln in parse_blast(handle):
        # extend block
        if (aln["qname"], aln["tname"]) == key:
            alns.append(aln)
        # yield block and start new one
        else:
            yield alns
            key = (aln["qname"], aln["tname"])
            alns = [aln]
    yield alns


def prune_alns(alns, min_length=0, min_evalue=1e-3):
    # remove alignments with < min_length or > min_evalue
    # discard alignments after the query length has been covered by 110%
    keep = []
    cur_aln = 0
    qry_len = alns[0]["qlen"]
    for aln in alns:
        qcoords = aln["qcoords"]
        aln_len = max(qcoords) - min(qcoords) + 1
        if aln_len < min_length or aln["evalue"] > min_evalue:
            continue
        if cur_aln >= qry_len or aln_len + cur_aln >= 1.10 * qry_len:
            break
        keep.append(aln)
        cur_aln += aln_len
    return keep


def compute_cov(alns):

    # merge qcoords
    coords = sorted([a["qcoords"] for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:

        # overlapping, update start coord
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)

        # non-overlapping, append to list
        else:
            nr_coords.append([start, stop])

    # compute query cov
    alen = sum([stop - start + 1 for start, stop in nr_coords])
    qcov = round(100.0 * alen / alns[0]["qlen"], 2)

    # merge tcoords
    coords = sorted([a["tcoords"] for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:
        # overlapping, update start coord
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)
        # non-overlapping, append to list
        else:
            nr_coords.append([start, stop])

    # recompute query cov
    alen = sum([stop - start + 1 for start, stop in nr_coords])
    tcov = round(100.0 * alen / alns[0]["tlen"], 2)

    return qcov, tcov


def ani_calculator(inpath, outpath):
    with open(outpath, "w") as out:
        fields = ["query", "reference", "ani", "qcov", "tcov", "norm_score"]
        out.write("\t".join(fields) + "\n")
        for alns in yield_alignment_blocks(open(inpath)):
            if alns is not None:
                alns = prune_alns(alns)
                if len(alns) > 0:
                    qname, tname = alns[0]["qname"], alns[0]["tname"]
                    ani = round(sum(a["len"] * a["pid"] for a in alns) / sum(
                        a["len"] for a in alns
                    ), 2)
                    qcov, tcov = compute_cov(alns)
                    norm_score = float(qcov) * ani / 100
                    row = [qname, tname, ani, qcov, tcov, norm_score]
                    out.write("\t".join([str(_) for _ in row]) + "\n")


##
## code to calculate AAI
##

def yield_diamond_hits(diamond):
    with open(diamond) as f:
        try :
            hits = [next(f).split()]
        except StopIteration:
            return
        for l in f:
            r = l.split()
            query = r[0].rsplit("_", 1)[0]
            last = hits[-1][0].rsplit("_", 1)[0]
            if query != last:
                yield last, hits
                hits = []
            hits.append(r)
        if len(hits) > 0:
            last = hits[-1][0].rsplit("_", 1)[0]
            yield last, hits


def split_hits(hits):
    target_to_hits = defaultdict(list)
    for hit in hits:
        tname = hit[1].rsplit("_", 1)[0]
        target_to_hits[tname].append(hit)
    return target_to_hits


def best_blast_hits(hits, query_key=0, score_key=-1):
    bhits = {}
    for hit in hits:
        if hit[query_key] not in bhits:
            bhits[hit[query_key]] = hit
        elif float(hit[score_key]) > float(bhits[hit[query_key]][score_key]):
            bhits[hit[query_key]] = hit
    return list(bhits.values())


def aai_main(inpath, outpath, selfpath):
    selfaai = {}
    for r in csv.DictReader(open(selfpath), delimiter="\t"):
        selfaai[r["genome_id"]] = float(r["selfscore"])
    with open(outpath, "w") as out:
        header = ["query", "reference", "hits", "aai", "raw_score", "norm_score"]
        out.write("\t".join(header) + "\n")
        for qname, hits in yield_diamond_hits(inpath):
            target_to_hits = split_hits(hits)
            for tname, thits in target_to_hits.items():
                bhits = best_blast_hits(thits)
                aai = mean([float(_[2]) for _ in bhits])
                score = sum([float(_[-1]) for _ in bhits])
                norm = 100 * score / selfaai[qname]
                row = [
                    qname,
                    tname,
                    len(bhits),
                    round(aai, 2),
                    round(score, 2),
                    round(norm, 2),
                ]
                out.write("\t".join([str(_) for _ in row]) + "\n")
