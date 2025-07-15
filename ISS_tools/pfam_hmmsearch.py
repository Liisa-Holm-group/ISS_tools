#!/usr/bin/env python3

import sys,os,re
import subprocess
import tempfile

import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from ISS_tools.pfam_io import parse_hmmer_output, write_FASTA
from ISS_tools.Pfam import Pfam

def add_pfam_annotation(df, seq_col, pfamdir, verbose=True):
    """
    Annotate input DataFrame with Pfam domains using a selected sequence column.

    Only the best Pfam hit per sbjct is retained!

    Args:
        df (pd.DataFrame): Input DataFrame with a sequence column (e.g., 'sequ-pileup').
        seq_col (str): Name of the column containing sequences to scan with hmmsearch.
        hmm_file (str): Path to Pfam-A.hmm database file.
        clan_file (str): Path to clan_pfam.tsv mapping file.
        verbose (bool): If True, prints progress messages to stderr.

    Returns:
        pd.DataFrame: Annotated DataFrame with Pfam information merged on 'sbjct'.
    """
    # use hard-coded hmm_file, clan_file
    pf=Pfam(pfamdir)
    # get hard-coded hmm_file and clan_file from Pfam module
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as tmp_fasta:
        fasta_path = tmp_fasta.name
        write_FASTA(df, seq_col=seq_col, outfile=fasta_path)

    try:
        if verbose:
            print(f"Running Pfam scan on sequences from column '{seq_col}'...", file=sys.stderr)

        parsed_df = run_pfam_scan(fasta_path, pf.hmm_file )
        # Step 3: Pick first hit per Sbjct
        parsed_df = parsed_df.sort_values(by='hmmer-evalue', ascending=True)
        filtered_df=parsed_df.drop_duplicates(subset=['sbjct'])
        clan_df = add_clan_annotation(filtered_df, pf.clan_file)
        result = pd.merge(df, clan_df, on='sbjct', how='left')

    finally:
        os.remove(fasta_path)

    return result

def add_clan_annotation(df, clan_file):
    # add CLAN annotation
    clan_df = pd.read_csv(clan_file, sep="\t")
    result = pd.merge(df, clan_df, left_on="pfam", right_on='AC', how="left")[['sbjct','pfam','CL','hmmer-evalue']]
    result.rename(columns={'CL':'clan'},inplace=True)
    # unassigned clan -> pfam
    result['clan']=result['clan'].fillna(result['pfam'])

    return(result)

def run_pfam_scan(fasta_file, hmm_file):
    """
    Run Pfam domain scan using hmmsearch and return results as a pandas DataFrame.
    Args:
        fasta_file (str): Path to input FASTA file.
        hmm_file (str): Path to Pfam-A.hmm database.
     Returns:
        pandas.DataFrame: DataFrame with columns ['sbjct', 'pfam', 'hmmmer_evalue', ...]
    """
    # Create temporary files
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmpout1, \
         tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmpout2:
        try:
            # Step 1: Run hmmsearch
            subprocess.run([
                "hmmsearch",
                "-o", tmpout1.name,
                "--acc", "--noali", "--notextw", "--cut_tc",
                hmm_file,
                fasta_file
            ], check=True)
            # Step 2: Parse HMMER output to DataFrame
            with open(tmpout1.name) as f:
               txt=f.readlines()
            parsed_df = parse_hmmer_output(txt)
        finally:
            os.remove(tmpout1.name)
            os.remove(tmpout2.name)
    return parsed_df

def pfam_scan(fasta_input, pfamdir, output_path, source_label=""):
    """
    Handle the full scan pipeline: prepares FASTA, runs hmmsearch, adds clan, and writes tabular output.

    Args:
        fasta_input: Either a Path to a FASTA file, or a pandas DataFrame with ['sbjct', 'sbjct-sequence'].
        pfamdir: Path to Pfam data directory.
        output_path: Output TSV file.
        source_label: String label for status messages.
    """
    pf=Pfam(pfamdir)
    #hmm_file = os.path.join(pfamdir, "Pfam-A.hmm")
    #clan_file = os.path.join(pfamdir, "Pfam_data.tsv")

    # Create a temporary FASTA file if input is a DataFrame
    if isinstance(fasta_input, pd.DataFrame):
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as tmp_fa:
            write_FASTA(fasta_input, outfile=tmp_fa.name)
            fasta_file = tmp_fa.name
            temp_created = True
    else:
        fasta_file = str(fasta_input)
        temp_created = False

    print(f"Running hmmsearch on {source_label}...")
    parsed_df = run_pfam_scan(fasta_file, pf.hmm_file)

    # add CLAN annotation
    result=add_clan_annotation(parsed_df, pf.clan_file)

    print(f"Writing final results to {output_path}")
    result.to_csv(output_path, sep="\t", index=False)

    if temp_created:
        os.remove(fasta_file)
