#!/usr/bin/env python3
import argparse
import sys
import os
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from ISS_tools.pfam_hmmsearch import add_pfam_annotation
from ISS_tools.parse_dali_txt import dali_txt_to_df, seriate


def main():
    parser = argparse.ArgumentParser(
        description="Post-process DALI output file (raw DALI .txt file or parsed .tsv)"
    )
    parser.add_argument(
        "--in", dest="infile", default=None, help="Input file (default: stdin)"
    )
    parser.add_argument(
        "--out", dest="outfile", default=None, help="Output file (default: stdout)"
    )
    parser.add_argument(
        "--parse", action="store_true", help="Parse DALI .txt file into TSV format"
    )
    parser.add_argument(
        "--annotate",
        action="store_true",
        help="Add Pfam annotations of aligned fragments",
    )
    parser.add_argument("--dat1", help="Path to DAT1 file (required with --parse)")
    parser.add_argument("--dat2", help="Path to DAT2 file (required with --parse)")
    parser.add_argument(
        "--target",
        dest="target",
        default="USER",
        help="Name of target database (optional with --parse)",
    )
    parser.add_argument("--seriate", action="store_true", help="Reorder rows prettily (slow on big data setes)")
    parser.add_argument(
        "--pfamdir", help="Path to HMM profile library (required with --annotate)"
    )

    args = parser.parse_args()

    # Print help if no arguments are given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Determine input source
    if args.infile:
        if args.parse and not args.infile.endswith(".txt"):
            sys.exit("Error: --in file must end with .txt when using --parse")
        elif not args.parse and not args.infile.endswith(".tsv"):
            sys.exit("Error: --in file must end with .tsv unless using --parse")
        infile = open(args.infile, "r")
    else:
        infile = sys.stdin

    # Parse if requested
    if args.parse:
        if not args.infile or not args.dat1 or not args.dat2:
            sys.exit("Error: --infile, --dat1 and --dat2 are required with --parse")
        df = dali_txt_to_df(args.infile, args.dat1, args.dat2, args.target)
    else:
        df = pd.read_csv(infile, sep="\t")

    if infile is not sys.stdin:
        infile.close()

    # clusternap order
    if args.seriate:
        df = seriate(df)

    # Pfam annotation if requested
    if args.annotate:
        if not args.pfamdir:
            sys.exit("Error: --pfamdir required with --annotate")
        df = add_pfam_annotation(df, "sequ-pileup", args.pfamdir)

    # Output
    if args.outfile:
        df.to_csv(args.outfile, sep="\t", index=False)
    else:
        df.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
