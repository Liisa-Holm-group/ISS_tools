#!/usr/bin/env python3

import sys,os
import subprocess
import tempfile
import urllib.request
from pathlib import Path
import gzip
import shutil 

import pandas as pd

### pfam_download ###
def parse_all_gf_blocks_from_gzip(gzip_path):
    """
    Parse all #=GF blocks from a gzipped Stockholm-format file into a DataFrame.

    Parameters:
        gzip_path (str): Path to a .gz file (e.g. Pfam-A.hmm.dat.gz)

    Returns:
        pd.DataFrame: One row per GF block; columns are GF field names.
    """
    with gzip.open(gzip_path, 'rt', encoding='utf-8') as f:
        content = f.read()

    blocks = content.strip().split('//')
    rows = []

    for block in blocks:
        gf_data = {}
        for line in block.strip().splitlines():
            if line.startswith("#=GF"):
                parts = line.split(maxsplit=2)
                if len(parts) == 3:
                    _, key, value = parts
                    gf_data[key] = value
        if gf_data:
            rows.append(gf_data)

    df=pd.DataFrame(rows)

    # strip version fro 'AC'
    df['AC_version'] = df['AC']
    df['AC'] = df['AC'].str.split('.').str[0]

    # rename
    #df.rename(columns={'AC': 'pfam', 'CL': 'clan'}, inplace=True)

    return df

def download_pfam(pfamdir):
    base_url = 'http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/'
    pfamdir = Path(pfamdir)
    pfamdir.mkdir(parents=True, exist_ok=True)

    files_to_download = ['Pfam-A.hmm.dat.gz', 'Pfam-A.hmm.gz']
    downloaded_files = {}

    for filename in files_to_download:
        url = base_url + filename
        local_gz_path = pfamdir / filename
        print(f"Downloading {filename} fron {url} to {local_gz_path}")
        urllib.request.urlretrieve(url, local_gz_path)
        downloaded_files[filename] = local_gz_path

    # Unzip Pfam-A.hmm.gz only
    hmm_gz = downloaded_files['Pfam-A.hmm.gz']
    hmm_path = hmm_gz.with_suffix('')  # Removes the .gz suffix
    print(f"Unzipping {hmm_gz} to {hmm_path}")
    with gzip.open(hmm_gz, 'rb') as f_in, open(hmm_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    # run hmmpress on downloaded hmm-library
    try:
            subprocess.run(["hmmpress", str(hmm_path)], check=True)
            print("hmmpress completed successfully.")
    except subprocess.CalledProcessError:
            print("Error: hmmpress failed. Make sure HMMER is installed and in your PATH.")

    # Parse the gzipped file
    df = parse_all_gf_blocks_from_gzip(os.path.join(pfamdir,"Pfam-A.hmm.dat.gz"))
    outfile = pfamdir / 'Pfam_data.tsv'
    df.to_csv(outfile, sep='\t', index=False)
    print(f"Wrote Pfam data to {outfile}")

    print("Download and extraction complete.")
    return downloaded_files


def ensure_hmmpress(hmm_file):
    """
    Run hmmpress on the given .hmm file if pressed index files are missing.
    """
    hmm_path = Path(hmm_file)
    expected_suffixes = [".h3f", ".h3i", ".h3m", ".h3p"]
    missing = any(not (hmm_path.with_suffix(suffix)).exists() for suffix in expected_suffixes)

    if missing:
        print(f"Running hmmpress on {hmm_file}...")
        try:
            subprocess.run(["hmmpress", str(hmm_file)], check=True)
            print("hmmpress completed successfully.")
        except subprocess.CalledProcessError:
            print("Error: hmmpress failed. Make sure HMMER is installed and in your PATH.")
    else:
        print("hmmpress index files already exist — skipping.")
