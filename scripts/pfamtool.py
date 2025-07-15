#!/usr/bin/env python3

import argparse
import sys, os
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from ISS_tools.pfam_download import download_pfam
from ISS_tools.pfam_io import dalidatdir_sequences_to_df, parse_hmmer_output
from ISS_tools.helpers import validate_required_args, check_program_exists
from ISS_tools.pfam_hmmsearch import pfam_scan

# action handler functions

def handle_download(args):
    validate_required_args(['pfamdir'], args)
    print("Downloading Pfam data...")
    return download_pfam(args.pfamdir)

def handle_search(args):
    validate_required_args(['pfamdir', 'search', 'out'], args)
    fasta_path = Path(args.search)
    pfam_scan(fasta_path, args.pfamdir, args.out, source_label="provided FASTA input")

def handle_parse(args):
    validate_required_args([], args)
    with open(args.parse,'r') as f:
        txt=f.readlines()
    df = parse_hmmer_output(txt)
    df.to_csv(args.out, sep='\t', index=False)

def handle_update(args):
    validate_required_args(['dat', 'pfamdir', 'out'], args)
    print("Extracting sequences from DALI .dat files...")
    seq_df = dalidatdir_sequences_to_df(args.dat)
    pfam_scan(seq_df, args.pfamdir, args.out, source_label="DALI-derived sequences")

def main():
    usage_examples = """\

Examples:

  Download hmm profile library of current Pfam release:
    python pfamtool.py --download --pfamdir ./pfamdata/

  Run hmmsearch on all sequences of proteins imported to DALI and stored in ./DAT:
    python pfamtool.py --update --dat ./DAT/ --pfamdir ./pfamdata/ --out ./pfadata/PDB.clan_hmmer_tc.tsv

  Run hmmsearch with a FASTA file and parse output:
    python pfamtool.py --search myseqs.fasta --pfamdir ./pfamdata/ --out hmmsearch_results.tsv

  Parse existing hmmsearch output:
    python pfamtool.py --parse results.txt --out parsed.tsv

"""
    
    parser = argparse.ArgumentParser(
        description="Pfam data utility script.\n\nhmmsearch is required with --update and --search.",
        epilog=usage_examples,
       formatter_class=argparse.RawTextHelpFormatter
    )

    parser = argparse.ArgumentParser(description="Pfam data utility script.\n\nhmmsearch is required with --update and --search.\n")
    parser.add_argument('--download', action='store_true', help='Download Pfam-A.hmm and Pfam-A.hmm.dat. Index and parse to TSV.')
    parser.add_argument('--update', action='store_true', help='Scan DALI .dat files and run hmmsearch on the sequences.')
    parser.add_argument('--pfamdir', type=str, default='../pfamdata/', help='Directory for Pfam data. Required with --download, --update, --search.')
    parser.add_argument('--dat', type=str, default='./DAT/', help='Directory containing DALI .dat files. Required with --update.')
    parser.add_argument('--out', type=str, default=sys.stdout, help='Output file for results. Defaults to STDOUT.')
    parser.add_argument('--search', type=str, help='Input FASTA file for hmmsearch. Parse results to TSV')
    parser.add_argument('--parse', type=str, help='Parse hmmsearch output file (txt) into TSV.')

    args = parser.parse_args()

    # Print help if no arguments provided
    if len(sys.argv) == 1:
        parser.print_help()
        # epilog not working...
        print(usage_examples)
        sys.exit(1)

    if args.download:
        args._action = 'download'
        handle_download(args)

    if args.search:
        check_program_exists('hmmsearch')
        args._action = 'search'
        handle_search(args)

    if args.parse:
        args._action = 'parse'
        handle_parse(args)

    if args.update:
        check_program_exists('hmmsearch')
        args._action = 'update'
        handle_update(args)

if __name__ == "__main__":
    main()


