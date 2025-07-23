import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pandas as pd


class Pfam:
    """
    Pfam database represented as pandas DataFrames.

    Attributes
    ----------
    hmm_file : Pfam-A profile library
    clan_file : file mapping Pfam families ('AC') to clans ('CL') and other information
    pfam_data_df : pandas.DataFrame
        DataFrame with Pfam family, clan, and descriptive information, including columns 'pfam', 'clan', and 'name'.

    Parameters
    ----------
    pfamdir : folder where Pfam was downloaded
    TRUE_df : Pandas DataFrame holding the mappings of target proteins to Pfam families and clans

    """

    def __init__(self, pfamdir):
        # Pfam download
        self.clan_file = f"{pfamdir}/Pfam_data.tsv"
        self.hmm_file = f"{pfamdir}/Pfam-A.hmm"
        # Pfam family information
        self.pfam_data_df = pd.read_csv(self.clan_file, sep="\t")
        # target database mapped to Pfam
        self.TRUE_df = None

    def load_TRUE(self, TRUE_TSV):
        """load target database to Pfam mappings
        Parameters:
               -----------
        TRUE_TSV : File with columns 'sbjct','pfam','clan'
        """
        self.TRUE_df = pd.read_csv(TRUE_TSV, sep="\t")

    def clan_member_families(self, CLAN):
        return self.pfam_data_df[self.pfam_data_df["CL"] == CLAN]

    def clan_member_proteins(self, CLAN):
        if self.TRUE_df is None:
            print(
                "WARNING: Pfam.TRUE_df is None - load it with Pfam.load_TRUE(TRUE_TSV)"
            )
            return pd.DataFrame
        else:
            return self.TRUE_df[self.TRUE_df["clan"] == CLAN]

    def add_TRUE(self, df, CLAN):
        "add TRUE column to df"
        if self.TRUE_df is None:
            print(
                "WARNING: Pfam.TRUE_df is None - load it with Pfam.load_TRUE(TRUE_TSV)"
            )
            df["TRUE"] = "unassigned"
        else:
            x_df = self.TRUE_df.loc[self.TRUE_df["clan"] == CLAN, ["clan", "sbjct"]]
            x_df = x_df.rename(columns={"clan": "TRUE"})
            new_df = pd.merge(df, x_df, on="sbjct", how="left")
            df["TRUE"] = new_df["TRUE"].astype(object).fillna("unassigned")
        return df

    def pfam_composition(self, df):
        """
        Annotate a DataFrame with Pfam and clan composition lists per subject.

        Parameters
        ----------
        df : pandas.DataFrame
            Input DataFrame containing at least a 'sbjct' column and related columns such as
            'z-score', 'ali-length', 'sbjct-length', 'description', and optionally 'motif'.

        Returns
        -------
        pandas.DataFrame
            A DataFrame merged with the original `df` containing two new columns:
            - 'pfamlist': list of unique Pfam entries per subject.
            - 'clanlist': list of unique clan entries per subject.

        Notes
        -----
        - The function groups Pfam and clan entries by 'sbjct', collecting unique values into lists.
        - Rows with no Pfam or clan assignments receive empty lists in the respective columns.
        """

        # Copy DataFrame
        sbjct_df = df[["sbjct"]].copy()

        # Left join with TRUE_df to get only relevant rows
        merged_df = pd.merge(
            sbjct_df, self.TRUE_df[["sbjct", "pfam", "clan"]], on="sbjct", how="left"
        )

        # Group and aggregate pfam and clan
        result_df = (
            merged_df.groupby("sbjct")
            .agg(
                {
                    "pfam": lambda x: list(set(x.dropna())) if x.notna().any() else [],
                    "clan": lambda x: list(set(x.dropna())) if x.notna().any() else [],
                }
            )
            .rename(columns={"pfam": "pfamlist", "clan": "clanlist"})
            .reset_index()
        )

        return pd.merge(result_df, df, on="sbjct", how="left")
