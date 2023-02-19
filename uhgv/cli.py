import argparse
import sys
import uhgv


def cli():
    parser = argparse.ArgumentParser(
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"""UHGV-toolkit v{uhgv.__version__}
https://github.com/snayfach/UHGV

usage: uhgv-tools <command> [options]

programs:
    download   download the UHGV genome database required for classify module
    classify   classify new genomes into phylogenetic groups from the UHGV""",
    )

    subparsers = parser.add_subparsers(help=argparse.SUPPRESS)

    download_database_parser = subparsers.add_parser(
        "download",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Download the database needed for classify module
\nusage: uhgv-tools download <destination>""",
    )
    uhgv.download.fetch_arguments(download_database_parser)

    classify_parser = subparsers.add_parser(
        "classify",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Classify new genomes into phylogenetic groups from the UHGV
\nusage: uhgv-tools classify [options]""",
    )
    uhgv.classify.fetch_arguments(classify_parser)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) == 2:
        if sys.argv[1] == "download":
            download_database_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "classify":
            classify_parser.print_help()
            sys.exit(0)

    args = vars(parser.parse_args())
    args["func"](args)
