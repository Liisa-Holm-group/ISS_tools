import sys,os,re,glob

import pandas as pd

def dalidatdir_sequences_to_df(dalidatdir) -> pd.DataFrame:
    records = []

    dat_files = glob.glob(os.path.join(dalidatdir, '*.dat'))

    for filepath in dat_files:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('-sequ'):
                    # Extract sequence after last quote
                    if '"' in line:
                        seq = line.strip().split('"')[-1]
                        # Get file stem and extract 5 characters before '.dat'
                        filename = os.path.basename(filepath)
                        dot_index = filename.rfind('.dat')
                        cd1 = filename[max(0, dot_index - 5):dot_index]
                        records.append({'sbjct': cd1, 'sbjct-sequence': seq})
                    break  # only first '-sequ' line needed

    return pd.DataFrame(records, columns=['sbjct', 'sbjct-sequence'])

def write_FASTA(df,id_col='sbjct',seq_col='sbjct-sequence',outfile='full.fasta'):
    "write sbjct, sequence in FASTA format"
    # write FASTA file
    df_clean = df.dropna()
    with open(outfile, "w") as f:
        for _, row in df_clean.iterrows():
            f.write(f">{row[id_col]}\n{row[seq_col]}\n")

def parse_hmmer_output(fileobj):
    """
    Parses multiple HMMER blocks from a single file-like input stream.
    Extracts Pfam accession, E-value, and Sequence from each.
    Args:
        fileobj: A file-like object (e.g., open file, StringIO)
    Returns:
        pandas.DataFrame with columns ['pfam', 'evalue', 'sbjct']
    """
    results = []
    pfam = None
    in_score_block = False
    skip_header_lines = 0
    for line in fileobj:
        line = line.strip()
        # Capture Pfam accession
        if line.startswith("Accession:"):
            match = re.search(r'(PF\d+)', line)
            if match:
                pfam = match.group(1)
        # Detect beginning of scores section
        elif line.startswith("Scores for complete sequences"):
            in_score_block = True
            skip_header_lines = 3  # skip 3 lines after this
        elif in_score_block and skip_header_lines > 0:
            skip_header_lines -= 1
        elif in_score_block:
            if not line:  # blank line: end of score block
                in_score_block = False
                continue
            # Parse a result line (ignore headers)
            tokens = line.split()
            if len(tokens) >= 9 and re.match(r'^[\deE.-]+$', tokens[0]):
                evalue = tokens[0]
                sequence = tokens[8]
                results.append((pfam, evalue, sequence))
        # End of one record
        elif line.startswith("//"):
            pfam = None
            in_score_block = False
            skip_header_lines = 0
    return pd.DataFrame(results, columns=["pfam", "hmmer-evalue", "sbjct"])

