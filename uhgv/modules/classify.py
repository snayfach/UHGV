#!/usr/bin/env python

import csv
import os
import sys
import math
import uhgv
import Bio.SeqIO
import shutil
import time
import gzip
import numpy as np
import subprocess as sp
from uhgv import utility
from collections import OrderedDict


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="classify")
    parser.add_argument(
        "-i",
        dest="input",
        type=str,
        required=True,
        metavar="PATH",
        help="Path to nucleotide seqs",
    )
    parser.add_argument(
        "-o",
        dest="outdir",
        type=str,
        required=True,
        metavar="PATH",
        help="Path to output directory",
    )
    parser.add_argument(
        "-d",
        dest="dbdir",
        type=str,
        required=True,
        metavar="PATH",
        help="Path to database directory",
    )
    parser.add_argument(
        "-s",
        dest="sens",
        choices=["fast", "sensitive", "very-sensitive"],
        default="sensitive",
        help="DIAMOND search sensitivity (sensitive)",
    )
    parser.add_argument(
        "-t",
        dest="threads",
        type=int,
        default=1,
        metavar="INT",
        help="Number of threads to run program with (1)",
    )
    parser.add_argument(
        "-p",
        dest="splits",
        type=int,
        metavar="INT",
        help="Number BLASTN jobs to spawn in parallel; 5GB RAM needed per job (default='threads')",
    )
    parser.add_argument(
        "--continue",
        action="store_true",
        default=False,
        help="Continue where program left off",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Suppress logging messages",
    )


def check_executables(requirements):
    fails = 0
    for program in requirements:
        found = False
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path.strip('"'), program)
            if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
                found = True
                break
        if not found:
            msg = f"\nError: required program '{program}' not executable or not found on $PATH"
            sys.stderr.write(msg)
            fails += 1
    if fails > 0:
        sys.exit("")


def blast_to_tophits(
    inpath, querykey="query", refkey="reference", scorekey="norm_score", refs=None
):
    """
    takes
        inpath: tsv file with header and at least 3 fields: query, reference, score
        querykey, refkey, scorekey: names of respective fields
        keep: optional list of references to keep, excluding the rest
    returns
        tophits: dict mapping each query to highest scoring record
    """
    tophits = {}
    for r in csv.DictReader(open(inpath), delimiter="\t"):
        r[scorekey] = float(r[scorekey])
        if refs and r[refkey] not in refs:
            continue
        elif r[querykey] not in tophits:
            tophits[r[querykey]] = r
        elif r[scorekey] > tophits[r[querykey]][scorekey]:
            tophits[r[querykey]] = r
    return tophits


def assign_aai_taxonomy(taxonomy, score):
    cutoffs = {"species": 95.0, "subgenus": 80.0, "genus": 65.0, "subfamily": 32.0, "family": 5.5}
    prefix2rank = {
        "vOTU": "species",
        "vSUBGEN": "subgenus",
        "vGENUS": "genus",
        "vSUBFAM": "subfamily",
        "vFAM": "family",
    }
    taxa = taxonomy.split(";")
    for i in range(2, len(taxa) + 1):
        lineage = taxa[: -i + 1]
        prefix = taxa[-i].split("-")[0]
        if prefix != "Unclassified":
            rank = prefix2rank[prefix]
            if score >= cutoffs[rank]:
                return ";".join(lineage)
    return None


class ViralClassifier:
    def __init__(self, args):
        self.set_args(args)
        self.set_paths()
        self.perform_checks()

    def set_args(self, args):
        self.args = args
        if self.args["splits"] is None:
            self.args["splits"] = self.args["threads"]
        self.args["blastcpus"] = math.ceil(1.0 * args["threads"] / args["splits"])

    def set_paths(self):
        self.paths = {}
        self.paths["input"] = self.args["input"]
        self.paths["outdir"] = self.args["outdir"]
        self.paths["dbdir"] = self.args["dbdir"].rstrip("/")
        self.paths["tmpdir"] = os.path.join(self.paths["outdir"], "tmp")
        self.paths["blastn"] = os.path.join(self.paths["tmpdir"], "blastn.tsv")
        self.paths["blastani"] = os.path.join(self.paths["tmpdir"], "blastani.tsv")
        self.paths["prodigal"] = os.path.join(self.paths["tmpdir"], "prodigal.faa")
        self.paths["selfaai"] = os.path.join(self.paths["tmpdir"], "selfscores.tsv")
        self.paths["diamond"] = os.path.join(self.paths["tmpdir"], "diamond.tsv")
        self.paths["blastaai"] = os.path.join(self.paths["tmpdir"], "blastaai.tsv")
        self.paths["blastdir"] = os.path.join(self.paths["tmpdir"], "blastn")
        self.paths["dmnddir"] = os.path.join(self.paths["tmpdir"], "diamond")
        self.paths["aaidir"] = os.path.join(self.paths["tmpdir"], "aai")
        self.paths["classify_summary"] = os.path.join(
            self.paths["outdir"], "classify_summary.tsv"
        )
        self.paths["taxon_info"] = os.path.join(
            self.paths["outdir"], "taxon_info.tsv"
        )

    def perform_checks(self):

        # check executables
        check_executables(["prodigal-gv", "diamond", "blastn"])

        # check database files
        if not os.path.exists(self.args["dbdir"]):
            sys.exit("\nError: database directory not found: '%s'" % args["dbdir"])
        files = [
            "genomes.fna",
            #"genomes.nal",
            "proteins.faa",
            "proteins.dmnd",
            "genome_taxonomy.tsv",
        ]
        for file in files:
            if not os.path.exists(os.path.join(self.paths["dbdir"], file)):
                sys.exit("\nError: database file not found: '%s'" % file)

        # check input seqs
        if not os.path.exists(self.args["input"]):
            sys.exit("\nError: input file does not exist")

        # check output directory does not exist
        # unless using --continue flag
        if os.path.exists(self.paths["tmpdir"]):
            if not self.args["continue"]:
                sys.exit(
                    "\nError: output directory already exists. Remove directory or use --continue"
                )
        else:
            os.makedirs(self.paths["tmpdir"])

    def load_queries(self):
        self.queries = OrderedDict()
        for r in Bio.SeqIO.parse(self.paths["input"], "fasta"):
            self.queries[r.id] = {}
            self.queries[r.id]["length"] = len(r.seq)

    def load_refdb(self):
        self.ref_genomes = {}
        path = os.path.join(self.paths["dbdir"], "genome_taxonomy.tsv")
        for r in csv.DictReader(open(path), delimiter="\t"):
            r["taxonomy"] = ";".join(
                [
                    r["family_vc"],
                    r["subfamily_vc"],
                    r["genus_vc"],
                    r["subgenus_vc"],
                    r["species_vc"],
                ]
            )
            r["taxonomy"] = r["taxonomy"].replace("NULL", "Unclassified")
            self.ref_genomes[r["genome_id"]] = r
        self.ref_clusters = {}
        path = os.path.join(self.paths["dbdir"], "viral_cluster_info.tsv")
        for r in csv.DictReader(open(path), delimiter="\t"):
            self.ref_clusters[r["taxon_id"]] = r
            del self.ref_clusters[r["taxon_id"]]["taxon_id"]

    def blastani(self):
        if os.path.exists(self.paths["blastani"]):
            return

        utility.split_fasta(
            self.paths["input"], self.paths["blastdir"], self.args["splits"], ".fna"
        )

        commands = []
        for file in os.listdir(self.paths["blastdir"]):
            if file.endswith(".fna"):
                inpath = os.path.join(self.paths["blastdir"], file)
                cmd = f"blastn -query {inpath} "
                cmd += f"-db {self.paths['dbdir']}/genomes "
                cmd += f"-out {self.paths['blastdir']}/{file}.tsv "
                cmd += f"-num_threads {self.args['blastcpus']} "
                cmd += "-outfmt '6 std qlen slen' "
                cmd += "-max_target_seqs 1000 "
                cmd += f"2> {self.paths['blastdir']}/{file}.log"
                commands.append([cmd])

        return_codes = utility.parallel(utility.run_shell, commands, self.args["splits"])
        if sum(return_codes) != 0:
            msg = "\nError: One or more blastn tasks failed to run\n"
            msg += f"See logs for details: {blastdir}/*.log"
            sys.exit(msg)

        with open(self.paths["blastn"], "w") as out:
            for file in os.listdir(self.paths["blastdir"]):
                if file.endswith(".tsv"):
                    inpath = os.path.join(self.paths["blastdir"], file)
                    for line in open(inpath):
                        out.write(line)

        utility.ani_calculator(self.paths["blastn"], self.paths["blastani"])

        shutil.rmtree(self.paths["blastdir"])

    def call_genes(self):
        if os.path.exists(self.paths["prodigal"]):
            return
        utility.parallel_prodigal(
            tmpdir=self.paths["tmpdir"],
            input=self.paths["input"],
            output=self.paths["prodigal"],
            threads=self.args["threads"],
            cleanup=True,
        )

    def self_protein_alignment(self):

        if not os.path.exists(self.paths["selfaai"]):

            data = {}

            # make dir
            selfaln_dir = os.path.join(self.paths["tmpdir"], "selfaln")
            if not os.path.exists(selfaln_dir):
                os.makedirs(selfaln_dir)

            # split fasta
            handle, qname, filenum = None, None, 0
            for r in Bio.SeqIO.parse(self.paths["prodigal"], "fasta"):
                genome = r.id.rsplit("_", 1)[0]
                if genome not in data:
                    data[genome] = {"genes": 0, "selfscore": 0}
                data[genome]["genes"] += 1
                if qname is None or genome != qname:
                    filenum += 1
                    qname = genome
                    handle = open(os.path.join(selfaln_dir, str(filenum)) + ".faa", "w")
                handle.write(">" + r.id + "\n" + str(r.seq) + "\n")
            handle.close()

            # run diamond
            commands = []
            files = [_ for _ in os.listdir(selfaln_dir) if _.endswith(".faa")]
            for file in files:
                path = os.path.join(selfaln_dir, file)
                cmd = "diamond blastp "
                cmd += "--masking none "
                cmd += "-k 1000 -e 1e-3 "
                cmd += f"--{self.args['sens']} "
                cmd += f"--query {path} "
                cmd += f"--db {path} "
                cmd += f"--out {path}.tsv "
                cmd += f"--threads 1 "
                cmd += f"2> {path}.log"
                commands.append([cmd])

            return_codes = utility.parallel(utility.run_shell, commands, self.args["splits"])
            if sum(return_codes) != 0:
                msg = "\nError: One or more diamond tasks failed to run\n"
                msg += f"See logs for details: {selfaln_dir}/*.log"
                sys.exit(msg)

            # calculate aai
            files = [_ for _ in os.listdir(selfaln_dir) if _.endswith(".tsv")]
            for file in files:
                path = os.path.join(selfaln_dir, file)
                score = 0
                genome = None
                for r in csv.reader(open(path), delimiter="\t"):
                    genome = r[0].rsplit("_", 1)[0]
                    if r[0] == r[1]:
                        score += float(r[-1])
                data[genome]["selfscore"] = score

            # write results
            with open(self.paths["selfaai"], "w") as out:
                header = ["genome_id", "genes", "selfscore"]
                out.write("\t".join(header) + "\n")
                for id in data:
                    rec = [id, data[id]["genes"], data[id]["selfscore"]]
                    out.write("\t".join([str(_) for _ in rec]) + "\n")

            # cleanup
            shutil.rmtree(selfaln_dir)

        # update queries
        for r in csv.DictReader(open(self.paths["selfaai"]), delimiter="\t"):
            self.queries[r["genome_id"]].update(r)

    def db_protein_alignment(self):
        if os.path.exists(self.paths["diamond"]):
            return
        cmd = "diamond blastp "
        cmd += "-k 1000 -e 1e-3 "
        cmd += "--masking none "
        cmd += f"--{self.args['sens']} "
        cmd += f"--query {self.paths['prodigal']} "
        cmd += f"--out {self.paths['diamond']} "
        cmd += f"--threads {self.args['threads']} "
        cmd += f"--db {self.paths['dbdir']}/proteins.dmnd "
        cmd += "--outfmt 6 "
        cmd += f"&> {self.paths['diamond']}.log"
        p = sp.Popen(cmd, shell=True)
        return_code = p.wait()
        if return_code != 0:
            msg = "\nError: DIAMOND database searched failed to run\n"
            msg += f"See log for details: {self.paths['diamond']}.log"
            sys.exit(msg)

    def blastaai(self):
        if os.path.exists(self.paths["blastaai"]):
            return
            
        utility.split_dmnd(self.paths["diamond"], self.paths["dmnddir"], self.args["threads"])
  
        if not os.path.exists(self.paths["aaidir"]):
            os.makedirs(self.paths["aaidir"])
             
        argument_list = []
        for file in os.listdir(self.paths["dmnddir"]):
            inpath = os.path.join(self.paths["dmnddir"], file)
            outpath = os.path.join(self.paths["aaidir"], file)
            argument_list.append([inpath, outpath, self.paths["selfaai"]])
        utility.parallel(utility.aai_main, argument_list, threads=self.args["threads"])

        with open(self.paths["blastaai"], "w") as out:
           
            for file in os.listdir(self.paths["aaidir"]):            
                handle = open(os.path.join(self.paths["aaidir"], file))
                out.write(next(handle))
                break
                        
            for file in os.listdir(self.paths["aaidir"]):            
                handle = open(os.path.join(self.paths["aaidir"], file))
                next(handle)
                for line in handle:
                    out.write(line)
                handle.close()
        
        shutil.rmtree(self.paths["dmnddir"])
        shutil.rmtree(self.paths["aaidir"])



    def find_top_hits(self):
        for type in ["blastani", "blastaai"]:
            tophits = blast_to_tophits(self.paths[type], refs=self.ref_genomes)
            for qname, hit in tophits.items():
                self.queries[qname][type] = hit

    def assign_taxonomy(self):
        for id in self.queries:

            r = {}
            r["genome_id"] = id
            r["genome_length"] = self.queries[id]["length"]
            ### temp fix
            r["genome_num_genes"] = self.queries[id]["genes"] if "genes" in self.queries[id] else "NA"
            r["taxon_id"] = None
            r["taxon_lineage"] = None
            r["class_method"] = None
            r["class_rank"] = None
            r["ani_reference"] = None
            r["ani_identity"] = None
            r["ani_query_af"] = None
            r["ani_target_af"] = None
            r["ani_taxonomy"] = None
            r["aai_reference"] = None
            r["aai_shared_genes"] = None
            r["aai_identity"] = None
            r["aai_score"] = None
            r["aai_taxonomy"] = None

            if "blastani" in self.queries[id]:
                h = self.queries[id]["blastani"]
                r["ani_reference"] = h["reference"]
                r["ani_identity"] = h["ani"]
                r["ani_query_af"] = h["qcov"]
                r["ani_target_af"] = h["tcov"]
                r["ani_taxonomy"] = self.ref_genomes[h["reference"]]["taxonomy"]

            if "blastaai" in self.queries[id]:
                h = self.queries[id]["blastaai"]
                r["aai_reference"] = h["reference"]
                r["aai_shared_genes"] = int(h["hits"])
                r["aai_identity"] = float(h["aai"])
                r["aai_score"] = float(h["norm_score"])
                r["aai_taxonomy"] = self.ref_genomes[h["reference"]]["taxonomy"]

            if (
                r["ani_reference"] is not None
                and float(r["ani_identity"]) >= 95
                and (float(r["ani_query_af"]) >= 85 or float(r["ani_target_af"]) >= 85)
            ):
                while r["ani_taxonomy"].endswith(";Unclassified"):
                    r["ani_taxonomy"] = r["ani_taxonomy"].rsplit(";Unclassified", 1)[0]
                r["taxon_lineage"] = r["ani_taxonomy"]
                r["taxon_id"] = r["taxon_lineage"].split(";")[-1]
                r["class_method"] = "ANI"
                r["class_rank"] = "species"

            elif r["aai_reference"] is not None:
                r["taxon_lineage"] = assign_aai_taxonomy(
                    r["aai_taxonomy"], r["aai_score"]
                )
                r["class_method"] = "AAI"
                if r["taxon_lineage"]:
                    rank_dict = {
                        "vFAM": "family",
                        "vSUBFAM": "subfamily",
                        "vGENUS": "genus",
                        "vSUBGEN": "subgenus",
                    }
                    r["taxon_id"] = r["taxon_lineage"].split(";")[-1]
                    r["class_rank"] = rank_dict[r["taxon_id"].split("-")[0]]

            if r["taxon_id"] is not None:
                r.update(self.ref_clusters[r["taxon_id"]])

            self.queries[id]["record"] = r

    def write_results(self):
        fields = [
            "genome_id",
            "genome_length",
            "genome_num_genes",
            "taxon_id",
            "class_method",
            "class_rank",
            "ani_reference",
            "ani_identity",
            "ani_query_af",
            "ani_target_af",
            "ani_taxonomy",
            "aai_reference",
            "aai_shared_genes",
            "aai_identity",
            "aai_score",
            "aai_taxonomy",
        ]
        with open(self.paths["classify_summary"], "w") as out:
            out.write("\t".join(fields) + "\n")
            for query in self.queries.values():
                rec = [query["record"][f] if f in query["record"] else "NA" for f in fields]
                out.write("\t".join([str(_) for _ in rec]) + "\n")

        fields = [
            "genome_id",
            "taxon_id",
            "taxon_lineage",
            "host_lineage",
            "ictv_lineage",
            "lifestyle",
            "genome_length_median",
            "genome_length_iqr",
#            "ncbi_genomes",
#            "ictv_genomes",
#            "uhgv_genomes",
#            "uhgv_complete",
#            "uhgv_high_quality",
#            "uhgv_medium_quality",
#            "uhgv_best_quality",
#            "ncbi_list",
        ]
        with open(self.paths["taxon_info"], "w") as out:
            out.write("\t".join(fields) + "\n")
            for query in self.queries.values():
                rec = [query["record"][f] if f in query["record"] else "NA" for f in fields]
                out.write("\t".join([str(_) for _ in rec]) + "\n")


######################
##      Main
######################


def main(args):

    prog_start = time.time()
    vclass = ViralClassifier(args)

    logger = utility.get_logger(args["quiet"])
    logger.info(f"\nUHGV-tools v{uhgv.__version__}: classify")

    logger.info("[1/10] Reading input sequences")
    vclass.load_queries()

    logger.info("[2/10] Reading database sequences")
    vclass.load_refdb()

    logger.info("[3/10] Estimating ANI with blastn")
    vclass.blastani()

    logger.info("[4/10] Identifying genes using prodigal-gv")
    vclass.call_genes()

    logger.info("[5/10] Performing self alignment")
    vclass.self_protein_alignment()

    logger.info("[6/10] Aligning proteins to database")
    vclass.db_protein_alignment()

    logger.info("[7/10] Calculating amino acid similarity scores")
    vclass.blastaai()

    logger.info("[8/10] Finding top database hits")
    vclass.find_top_hits()

    logger.info("[9/10] Performing phylogenetic assignment")
    vclass.assign_taxonomy()

    logger.info("[10/10] Writing output file(s)")
    vclass.write_results()

    logger.info("\nSuccess!")
    logger.info("Elapsed time (s): %s" % round(time.time() - prog_start, 2))
    logger.info("Peak RAM usage (GB): %s" % round(utility.max_mem_usage(), 2))
