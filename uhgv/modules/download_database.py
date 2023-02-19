import os
import sys
import shutil
import time
import urllib.request
import subprocess as sp
from uhgv import utility
import uhgv


class DatabaseDownloader:
    def __init__(self, destination):
        self.url = "https://portal.nersc.gov/UHGV/toolkit/"
        self.destination = destination
        self.version = (
            urllib.request.urlopen(self.url + "CURRENT_RELEASE.txt")
            .read()
            .decode("utf-8")
            .strip()
        )
        self.filename = self.version + ".tar.gz"
        self.output_file = os.path.join(self.destination, self.filename)

    def download(self):
        database_url = self.url + self.filename
        with urllib.request.urlopen(database_url) as response:
            with open(self.output_file, "wb") as fout:
                shutil.copyfileobj(response, fout)

    def extract(self):
        shutil.unpack_archive(self.output_file, self.destination, "gztar")
        os.remove(self.output_file)

    def blastn_makedb(self):
        self.dbdir = os.path.join(
            self.destination, self.filename.replace(".tar.gz", "")
        )
        cmd = "makeblastdb "
        cmd += f"-in {self.dbdir}/genomes.fna "
        cmd += f"-out {self.dbdir}/genomes "
        cmd += "-dbtype nucl "
        cmd += "1> /dev/null "
        cmd += f"2> {self.dbdir}/genomes.log "
        p = sp.Popen(cmd, shell=True)
        return_code = p.wait()
        if return_code != 0:
            msg = "\nError: BLASTN database failed to build\n"
            msg += f"See log for details: {self.dbdir}/genomes.log"
            sys.exit(msg)

    def diamond_makedb(self):
        self.dbdir = os.path.join(
            self.destination, self.filename.replace(".tar.gz", "")
        )
        cmd = "diamond makedb "
        cmd += f"--in {self.dbdir}/proteins.faa "
        cmd += f"--db {self.dbdir}/proteins "
        cmd += f"--threads 1 "
        cmd += "1> /dev/null "
        cmd += f"2> {self.dbdir}/proteins.log "
        p = sp.Popen(cmd, shell=True)
        return_code = p.wait()
        if return_code != 0:
            msg = "\nError: DIAMOND database failed to build\n"
            msg += f"See log for details: {self.dbdir}/proteins.log"
            sys.exit(msg)


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="download_database")
    parser.add_argument(
        "destination",
        type=str,
        help="Directory where the database will be downloaded to.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Suppress logging messages",
    )


def main(args):

    program_start = time.time()
    logger = utility.get_logger(args["quiet"])
    if not os.path.exists(args["destination"]):
        os.makedirs(args["destination"])

    logger.info(f"\nUHGV-tools v{uhgv.__version__}: download_database")

    logger.info("[1/5] Checking latest version of database...")
    db = DatabaseDownloader(args["destination"])

    logger.info(f"[2/5] Downloading '{db.version}'...")
    db.download()

    logger.info(f"[3/5] Extracting '{db.version}'...")
    db.extract()

    logger.info(f"[4/5] Building BLASTN database...")
    db.blastn_makedb()

    logger.info(f"[5/5] Building DIAMOND database...")
    db.diamond_makedb()

    logger.info("Run time: %s seconds" % round(time.time() - program_start, 2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(), 2))
