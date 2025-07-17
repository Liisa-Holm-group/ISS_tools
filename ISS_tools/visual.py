#!/usr/bin/env python
# coding: utf-8

# # Prepare environment

import os
import sys
import itertools


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from ISS_tools.helpers import import_install
from ISS_tools.parse_dali_txt import pileup

# Third-party packages
import_install("scipy")
import_install("numpy")
import_install("matplotlib")
import_install("pandas")
import_install("seaborn")
import_install("logomaker")
import_install("PIL")
import_install("bitstring")

from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns
import logomaker as lm

import_install("plotly")
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "notebook"  # or 'inline' or 'notebook_connected'

# # Setup Query


def read_tsv(tsvfile):
    "returns DALIFrame"
    df = DALIFrame()
    df.load_tsv(tsvfile)
    return df


def read_dali_tsv(tsvfile):
    "reads header-less TSV file assumed to be output of ISS_rundali.py"
    dali_tsv_colnames = [
        "cd1",
        "sbjct",
        "z-score",
        "rmsd",
        "ali-length",
        "sbjct-length-cropped",
        "seq-identity",
        "qstarts",
        "sstarts",
        "lengths",
        "rotation",
        "translation",
        "null1",
        "null2",
    ]
    df = pd.read_csv(tsvfile, sep="\t", header=None, names=dali_tsv_colnames)
    print(f"Read {df.shape[0]} rows with {df.shape[1]} columns from {tsvfile}")
    return df


def cleaning(df):
    # lowercase cols
    df.columns = [col.lower() for col in df.columns]
    # integers
    for col in ("ali-length", "sbjct-length"):
        df[col] = (
            pd.to_numeric(df[col], errors="coerce").astype("Int64").fillna(0)
        )  # nullable integer type
    # create pileups if issing
    nres = df.loc[df["database"] == "Query", "sbjct-length"].values[0]
    if "sbjct-coverage" not in df.columns:
        df["query-coverage"] = df["ali-length"] / nres
        df["sbjct-coverage"] = df["ali-length"] / nres
    # create coverages if issing
    if "dssp-pileup" not in df.columns:
        df["dssp-pileup"] = df.apply(
            lambda row: pileup(
                row["qstarts"], row["sstarts"], row["lengths"], row["sbjct-dssp"], nres
            ),
            axis=1,
        )
        df["sequ-pileup"] = df.apply(
            lambda row: pileup(
                row["qstarts"],
                row["sstarts"],
                row["lengths"],
                row["sbjct-sequence"],
                nres,
            ),
            axis=1,
        )
    # floats
    for col in ("z-score", "rmsd", "query-coverage", "sbjct-coverage", "seq-identity"):
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    return df


def add_TRUE(df, CLAN, TRUE_df):
    """Add a TRUE column indicating whether the sbjct belongs to the specified CLAN."""
    # Identify sbjct entries for the given CLAN
    match_sbjcts = TRUE_df.loc[TRUE_df["clan"] == CLAN, "sbjct"].unique()

    # Initialize TRUE column if it doesn't exist
    df["TRUE"] = df["sbjct"].apply(
        lambda x: CLAN if x in match_sbjcts else "unassigned"
    )

    return DALIFrame(df)


class DALIFrame:
    """
    Wrapper class around a pandas DataFrame with additional methods for handling DALI output.
    """

    def __init__(self, df=None):
        if df is not None:
            self._df = df.copy()
        else:
            self._df = pd.DataFrame()  # default to empty

    @classmethod
    def from_df(cls, df):
        """Wrap an existing DataFrame in a DALIFrame"""
        return cls(df)

    def add_TRUE(self, CLAN, TRUE_df):
        """Add a TRUE column indicating whether the sbjct belongs to the specified CLAN."""
        # Identify sbjct entries for the given CLAN
        match_sbjcts = TRUE_df.loc[TRUE_df["clan"] == CLAN, "sbjct"].unique()

        # Initialize TRUE column if it doesn't exist
        self._df["TRUE"] = self._df["sbjct"].apply(
            lambda x: CLAN if x in match_sbjcts else "unassigned"
        )

        # return self

    def subset(self, mask):
        """Return a DALIFrame with rows matching the mask.

        subset = df.subset((df['ali-length'] > 80) & (df['z-score'] > 0))
        """
        return DALIFrame(self._df[mask].copy())

    def sort_values(self, by, ascending=True, **kwargs):
        sorted_df = self._df.sort_values(by=by, ascending=ascending, **kwargs)
        return DALIFrame(sorted_df)

    def load_tsv(self, tsvfile):
        """
        Load TSV file and preprocess it.
        """
        self._df = cleaning(pd.read_csv(tsvfile, sep="\t"))
        self._standardize_columns()
        return self  # allow chaining

    def grouplabel(self, k=10, labelin="clan", labelout="hue"):
        """sort df by z-score, use top k labels, unassigned and pool others"""
        unique_clans = (
            self._df.sort_values("z-score", ascending=False)[labelin]
            .dropna()
            .drop_duplicates()
        )
        self._df[labelout] = "other"
        for x in unique_clans[:k]:
            print("assign", x)
            self._df.loc[self._df[labelin] == x, labelout] = self._df[labelin]

    def _standardize_columns(self):
        """
        Internal method to lowercase column names and clean up the DataFrame.
        """
        # Columns to convert and their desired types
        columns_to_convert = {
            "z-score": float,
            "rmsd": float,
            "query-coverage": float,
            "sbjct-coverage": float,
            "sbjct-length": int,
            "ali-length": int,
            "seq-identity": int,
        }

        # Convert each column to numeric, replacing '.' with 0 before casting
        for col, col_type in columns_to_convert.items():
            self._df[col] = self._df[col].replace(".", 0).astype(col_type)

        # Replace zero alignment lengths with the minimum non-zero length
        min_ali_len = self._df[self._df["ali-length"] > 0]["ali-length"].min()
        self._df["ali-length"] = self._df["ali-length"].replace(0, min_ali_len)

        # Fill missing Pfam and Clan annotations
        self._df["pfam"] = self._df["pfam"].astype(object).fillna("unassigned")
        self._df["clan"] = self._df["clan"].astype(object).fillna("unassigned")

        # Generate hovertext depending on available columns
        if "accession" in self._df.columns and "description" in self._df.columns:
            self._df["hovertext"] = (
                self._df["sbjct"]
                + "|"
                + self._df["accession"]
                + " "
                + self._df["description"]
            )
        else:
            self._df["hovertext"] = self._df["sbjct"] + " " + self._df["description"]

        # drop QUery row
        self._df.drop(self._df[self._df["database"] == "Query"].index, inplace=True)
        # self._df.drop(self._df[self._df['z-score'] < 2 ].index, inplace=True)

    def set_motif(self, positions, k=4, mute=False):
        """create motif column concatenating sequ-pileup letters in positions list
        print top-10 groups
        pool groups >k in group 'other'
        """
        # create motif column
        self._df["motif"] = self._df["sequ-pileup"].apply(
            lambda s: "".join([s[i] if len(s) > i else "-" for i in positions])
        )
        if not mute:
            print("\nFrequent")
            print(self._df["motif"].value_counts()[:10])
        # find top-k groups
        topk = self._df["motif"].value_counts().nlargest(k).index
        # print(topk)
        # reduce number of distinct groups
        self._df["motif"] = self._df["motif"].where(
            self._df["motif"].isin(topk), other="other"
        )
        if not mute:
            print("\nSelected")
            print(self._df["motif"].value_counts())

    def __getattr__(self, attr):
        """
        Delegate attribute access to the internal DataFrame.
        This allows calling DataFrame methods like df.head(), df['col'], etc.
        """
        return getattr(self._df, attr)

    def __getitem__(self, key):
        return self._df[key]

    def __setitem__(self, key, value):
        self._df[key] = value

    def __repr__(self):
        return repr(self._df)

    def __iter__(self):
        return iter(self._df)

    def scatter(
        self,
        x="z-score",
        y="ali-length",
        color=None,
        marker=None,
        left=None,
        right=None,
        bottom=None,
        top=None,
        title="Scatter Plot",
        hovertext=None,
        size=None,
    ):
        """interface to matplotlib scatterplots with continuous or categorical color,
        categorical markers, colorbar and marker legends
        """
        _scatter(
            self._df,
            x=x,
            y=y,
            color=color,
            marker=marker,
            left=left,
            right=right,
            bottom=bottom,
            top=top,
            title=title,
            size=size,
        )

    def iscatter(
        self,
        x="z-score",
        y="ali-length",
        color="seq-identity",
        marker=None,
        left=None,
        right=None,
        bottom=None,
        top=None,
        hovertext="hovertext",
        title=None,
        size=None,
    ):
        """interface to Plotly scatterplots with continuous or categorical color,
        categorical markers, colorbar and marker legends.
        """

        if color is None:
            print("Error: color must be defined")
            return
        else:
            series = self._df[color]
            color_is_cat = _is_categorical(series)
        if marker is not None:
            series = self._df[marker]
        if hovertext is not None:
            if hovertext not in self._df.columns:
                hovertext = None

        if not color_is_cat and marker is None:
            return _scatter_plotly(
                self._df,
                x=x,
                y=y,
                color=color,
                hovertext=hovertext,
                title=title,
                size=size,
            )
        if not color_is_cat and marker:
            return _scatter_plotly(
                self._df,
                x=x,
                y=y,
                color=color,
                marker=marker,
                hovertext=hovertext,
                title=title,
                size=size,
            )
        if color_is_cat and marker is None:
            return _dual_marker_plotly(
                self._df,
                x=x,
                y=y,
                color=color,
                marker=color,
                hovertext=hovertext,
                title=title,
                size=size,
            )
        if color_is_cat and marker:
            return _dual_marker_plotly(
                self._df,
                x=x,
                y=y,
                color=color,
                marker=marker,
                hovertext=hovertext,
                title=title,
                size=size,
            )
        return


def _is_categorical(series):
    """Infer if a pandas Series is categorical."""
    return series.dtype == "object" or isinstance(series, pd.CategoricalDtype)


def _scatter(
    df,
    x="z-score",
    y="ali-length",
    color=None,
    marker=None,
    left=None,
    right=None,
    bottom=None,
    top=None,
    title="Smart Scatter Plot",
    size=None,
):
    """
    Create a scatter plot from a DataFrame with flexible options for coloring, markers, and axis limits.

    Parameters:
    -----------
    df : pandas.DataFrame
        The DataFrame containing the data to plot.

    x : str, default='z-score'
        Column name in `df` to use for the x-axis values.

    y : str, default='ali-length'
        Column name in `df` to use for the y-axis values.

    color : str or None, optional
        Column name in `df` to determine the color of points. If None, all points have the same color.

    marker : str or None, optional
        Column name in `df` to determine the marker style of points. If None, all points use the default marker.

    left : float or None, optional
        Left boundary limit for the x-axis. If None, matplotlib chooses automatically.

    right : float or None, optional
        Right boundary limit for the x-axis. If None, matplotlib chooses automatically.

    bottom : float or None, optional
        Bottom boundary limit for the y-axis. If None, matplotlib chooses automatically.

    top : float or None, optional
        Top boundary limit for the y-axis. If None, matplotlib chooses automatically.

    size : str or None, optional
        Column name in df to determine the size of poitns. If None, all points have the same size.

    title : str, default='Smart Scatter Plot'
        Title of the plot.

    Returns:
    --------
    None

    Notes:
    ------
    - If `color` or `marker` correspond to categorical columns, the function will automatically assign colors and markers.
    - Axis limits can be set individually using `left`, `right`, `bottom`, and `top`.
    - This function simplifies creating informative scatter plots directly from DataFrame columns.
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.data = []  # tabula rasa
    default_color = "gray"
    default_marker = "o"

    # Marker style cycling
    # base_marker_styles = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'h', 'H', 'd', '|', '_']
    base_marker_styles = [
        "o",
        "s",
        "^",
        "D",
        "v",
        "P",
        "*",
        "X",
        "<",
        ">",
        "H",
        "h",
        "d",
    ]

    color_is_categorical = False
    marker_is_categorical = False

    # Handle color=None
    if color is None:
        if marker is None:
            ax.scatter(
                df[x],
                df[y],
                color=default_color,
                marker=default_marker,
                edgecolor="k",
                s=_point_size(df, size),
                alpha=0.7,
            )
        else:
            marker_vals = df[marker].unique()
            marker_cycle = itertools.cycle(base_marker_styles)
            marker_map = {m: next(marker_cycle) for m in marker_vals}

            for m in marker_vals:
                group = df[df[marker] == m]
                ax.scatter(
                    group[x],
                    group[y],
                    color=default_color,
                    marker=marker_map[m],
                    edgecolor="k",
                    s=_point_size(group, size),
                    alpha=0.7,
                )

            # Marker legend
            marker_legend = [
                Line2D(
                    [0],
                    [0],
                    marker=marker_map[m],
                    color="k",
                    label=m,
                    markerfacecolor="lightgray",
                    markersize=8,
                    linestyle="",
                )
                for m in marker_vals
            ]
            ax.legend(handles=marker_legend, title=marker, loc="lower right")

    # Handle color as continuous
    elif pd.api.types.is_numeric_dtype(df[color]):
        norm = plt.Normalize(vmin=df[color].min(), vmax=df[color].max())

        if marker is None:
            scatter = ax.scatter(
                df[x],
                df[y],
                c=df[color],
                cmap="viridis",
                marker=default_marker,
                edgecolor="k",
                s=_point_size(df, size),
                alpha=0.7,
            )
        else:
            marker_vals = df[marker].unique()
            marker_cycle = itertools.cycle(base_marker_styles)
            marker_map = {m: next(marker_cycle) for m in marker_vals}

            for m in marker_vals:
                group = df[df[marker] == m]
                ax.scatter(
                    group[x],
                    group[y],
                    c=group[color],
                    cmap="viridis",
                    norm=norm,
                    marker=marker_map[m],
                    edgecolor="k",
                    s=_point_size(group, size),
                    alpha=0.7,
                )

            # Marker legend
            marker_legend = [
                Line2D(
                    [0],
                    [0],
                    marker=marker_map[m],
                    color="k",
                    label=m,
                    markerfacecolor="gray",
                    markersize=8,
                    linestyle="",
                )
                for m in marker_vals
            ]
            ax.legend(handles=marker_legend, title=marker, loc="lower right")

        # Colorbar
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        plt.colorbar(
            sm, ax=ax, label=color, orientation="vertical", pad=0.02, location="right"
        )

    # Handle color as categorical
    else:
        color_is_categorical = True
        color_vals = df[color].unique()
        palette = sns.color_palette("tab10", len(color_vals))
        color_map = dict(zip(color_vals, palette))

        if marker is None:
            for c in color_vals:
                group = df[df[color] == c]
                ax.scatter(
                    group[x],
                    group[y],
                    color=color_map[c],
                    marker=default_marker,
                    edgecolor="k",
                    s=_point_size(group, size),
                    alpha=0.7,
                )
        else:
            marker_is_categorical = True
            marker_vals = df[marker].unique()
            marker_cycle = itertools.cycle(base_marker_styles)
            marker_map = {m: next(marker_cycle) for m in marker_vals}
            for m in marker_vals:
                for c in color_vals:
                    group = df[(df[marker] == m) & (df[color] == c)]
                    ax.scatter(
                        group[x],
                        group[y],
                        color=color_map[c],
                        marker=marker_map[m],
                        edgecolor="k",
                        s=_point_size(group, size),
                        alpha=0.7,
                    )

    # Custom legend for color (categorical)
    if color_is_categorical:
        color_legend = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label=cat,
                markerfacecolor=color_map[cat],
                markersize=8,
            )
            for cat in color_vals
        ]
        legend1 = ax.legend(handles=color_legend, title=color, loc="upper left")
        ax.add_artist(legend1)

    # Custom legend for marker (categorical)
    if marker_is_categorical:
        marker_legend = [
            Line2D(
                [0],
                [0],
                marker=marker_map[cat],
                color="k",
                label=cat,
                markerfacecolor="gray",
                markersize=8,
                linestyle="",
            )
            for cat in marker_vals
        ]
        legend2 = ax.legend(handles=marker_legend, title=marker, loc="lower right")

    # Optionally set horizontal plot limits
    if left is not None:
        plt.xlim(left=left)
    if right is not None:
        plt.xlim(right=right)
    if bottom is not None:
        plt.ylim(bottom=bottom)
    if top is not None:
        plt.ylim(top=top)

    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(title)
    plt.tight_layout()
    plt.show()
    plt.close()


# ## Scatterplot with interactive zooming and text on hover


def _scatter_plotly(
    df: pd.DataFrame,
    x: str = "z-score",
    y: str = "ali-length",
    color: str = None,
    hovertext: str = "hovertext",
    marker: str = None,  # Optional symbol-style grouping
    title: str = None,
    size: str = None,
) -> go.Figure:
    """
    Create an interactive scatter plot using Plotly from a DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame containing the data to plot.

    x : str, default='z-score'
        Column name to use for the x-axis values.

    y : str, default='ali-length'
        Column name to use for the y-axis values.

    color : str or None, optional
        Column name in `df` to map point colors. If None, all points use the same default color.

    hovertext : str, default='hovertext'
        Column name in `df` to display as hover text for each point.

    marker : str or None, optional
        Column name in `df` to determine marker symbols for grouping. If None, a single marker style is used.

    title : str or None, optional
        Title of the plot. If None, no title is set.

    Returns
    -------
    None

    Notes
    -----
    - Supports categorical coloring and symbol grouping.
    - Hover text enhances interactivity by showing additional information on mouse-over.
    """

    # If marker grouping is provided, we must split into multiple traces
    if marker:
        fig = go.Figure()

        # Define a list of symbols to cycle through
        marker_symbols = [
            "circle",
            "square",
            "diamond",
            "cross",
            "x",
            "triangle-up",
            "triangle-down",
        ]
        unique_markers = df[marker].unique()
        symbol_map = {
            m: marker_symbols[i % len(marker_symbols)]
            for i, m in enumerate(unique_markers)
        }

        for m in unique_markers:
            sub_df = df[df[marker] == m]
            trace_kwargs = dict(
                x=sub_df[x],
                y=sub_df[y],
                mode="markers",
                name=str(m),
                marker=dict(
                    color=sub_df[color],
                    colorscale="Viridis",
                    symbol=symbol_map[m],
                    size=_point_size(sub_df, size, sizedef="diameter"),
                    line_width=1,
                    showscale=False,  # We'll show the scale once outside the loop
                ),
            )
            if hovertext is not None:
                trace_kwargs["hovertext"] = sub_df[hovertext]
            fig.add_trace(go.Scattergl(**trace_kwargs))

        # Add a dummy colorbar using the full dataset
        fig.add_trace(
            go.Scattergl(
                x=[None],
                y=[None],  # invisible
                mode="markers",
                marker=dict(
                    color=df[color],
                    colorscale="Viridis",
                    colorbar=dict(title=color),
                    showscale=True,
                    size=0,  # invisible
                ),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    else:
        # Single-trace version
        trace = go.Scattergl(
            x=df[x],
            y=df[y],
            mode="markers",
            hovertext=df[hovertext],
            marker=dict(
                color=df[color],
                colorscale="Viridis",
                colorbar=dict(title=color),
                line_width=1,
                size=_point_size(df, size, sizedef="diameter"),
                showscale=True,
            ),
        )
        fig = go.Figure(trace)

    # Layout
    fig.update_layout(
        width=700,
        height=600,
        xaxis_title=x,
        yaxis_title=y,
        title=dict(text="Interactive ScatterGL Plot", x=0.5, xanchor="center"),
        template="plotly_white",
        legend=dict(
            x=1,
            y=0,
            xanchor="right",
            yanchor="bottom",
            bgcolor="rgba(255,255,255,0.6)",
            bordercolor="gray",
            borderwidth=1,
        ),
    )
    fig.show()


# # msa-heatmap+logo
class msa:
    """
    Visual representation of multiple sequence alignments of amino acid sequence or secondary structure.

    This class takes a DataFrame containing a column of aligned sequences or structural annotations
    (e.g., DSSP secondary structure codes or amino acid sequences) and creates a color-coded
    heatmap representation for visualization.

    Parameters
    ----------
    df : pandas.DataFrame

    col : str. Must be one of:
        - 'dssp-pileup' for secondary structure annotations (H, E, L)
        - 'sequ-pileup' for aligned amino acid sequences

    Notes
    -----
    Color mapping is hard-coded:
    - For 'dssp': H=blue, E=red, L=green
    - For 'aa': groups amino acids by physico-chemical properties with distinct colors
    """

    def __init__(self, df, col):
        self.df = df[df[col] != "."].dropna(subset=[col])
        if col == "dssp-pileup":
            self.pileup_col = "dssp-pileup"
            self.char_colors = self._dssp_colors()
            self.color_mapping = self._dssp_color_mapping()
        elif col == "sequ-pileup":
            self.pileup_col = "sequ-pileup"
            self.char_colors = self._aa_colors()
            self.color_mapping = self._aa_color_mapping()
        else:
            print("Error: col must be 'sequ-pileup' or 'dssp-pileup'")
        # initialize data
        self._rgb_array = None
        self._count_df = None
        # clan/pfam occupancy profiles
        self._grouped = None
        self._corr_matrix = None

    def _rgb(self):
        # Convert to RGB array
        matrix = self.df[self.pileup_col].astype(str).tolist()
        lengths = [len(row) for row in matrix]
        assert len(set(lengths)) == 1, f"Inconsistent row lengths: {lengths}"

        self._rgb_array = np.zeros(
            (self.df.shape[0], self.df.shape[1], 3), dtype=np.uint8
        )
        self._rgb_array = np.array(
            [
                [self.char_colors.get(char, (255, 255, 255)) for char in row]
                for row in matrix
            ]
        )

    def logo(self, logotype="information", left=None, right=None):
        "Displays the sequence logo using logomaker with optional axis limits."
        # Create a new figure with a fixed size
        plt.figure(figsize=(10, 3))

        # Generate logo matrix from the input dataframe and column
        if self._count_df is None:
            # Filter out invalid entries (e.g., single-character lines)
            seqs = (
                self.df[self.df[self.pileup_col].astype(str).str.len() > 1][
                    self.pileup_col
                ]
                .astype(str)
                .tolist()
            )
            subseqs = [s for s in seqs]
            # Convert sequences to a logomaker matrix (e.g., information content matrix)
            self.counts_df = lm.alignment_to_matrix(
                sequences=subseqs,
                to_type=logotype,
                characters_to_ignore=".-",  # Common gaps or unknowns in sequence alignments
            )

        # Create the sequence logo using Logomaker
        lm.Logo(self.counts_df, color_scheme=self.color_mapping)

        # Optionally set horizontal plot limits
        if left is not None:
            plt.xlim(left=left)
        if right is not None:
            plt.xlim(right=right)

        # Display the logo
        plt.show()

        # Close the plot to avoid reuse or memory leaks in Jupyter
        plt.close()

    def heatmap(
        self,
        xlabel="Position in sequence",
        ylabel="Row index",
        left=None,
        right=None,
        bottom=None,
        top=None,
        title=None,
    ):
        "Displays the RGB-based heatmap using matplotlib with optional axis labels and limits."
        if self._rgb_array is None:
            self._rgb()
        plt.figure(figsize=(6, 4), dpi=100)
        plt.imshow(self._rgb_array, aspect="auto")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if left is not None:
            plt.xlim(left=left)
        if right is not None:
            plt.xlim(right=right)
        if title is not None:
            plt.title(title)
        if top is not None:
            plt.ylim(top=top)
        if bottom is not None:
            plt.ylim(bottom=bottom)
        # Display the logo
        plt.show()
        # Close the plot to avoid reuse or memory leaks in Jupyter
        plt.close()

    def _dssp_colors(self):
        return {
            "H": (0, 0, 255),  # blue
            "E": (255, 0, 0),  # red
            "L": (0, 255, 0),  # green
        }

    def _aa_colors(self):
        return {
            # Hydrophobic - green/light blue
            "A": (102, 204, 102),  # soft green
            "V": (102, 204, 102),
            "L": (102, 204, 102),
            "I": (102, 204, 102),
            "M": (102, 204, 102),
            "P": (153, 204, 255),  # light blue
            "G": (153, 204, 255),
            # Cysteine - yellow
            "C": (255, 255, 102),  # bright yellow
            # Negatively charged (acidic) - red
            "D": (255, 51, 51),
            "E": (255, 51, 51),
            # Positively charged (basic) - blue
            "K": (51, 102, 255),
            "R": (51, 102, 255),
            "H": (102, 153, 255),
            # Aromatic - magenta
            "F": (255, 102, 255),
            "Y": (255, 102, 255),
            "W": (255, 102, 255),
            # Polar uncharged (green)
            "S": (0, 200, 0),
            "T": (0, 200, 0),
            "N": (0, 200, 0),
            "Q": (0, 200, 0),
        }

    def _dssp_color_mapping(self):
        return {"H": "blue", "E": "red", "L": "green"}

    def _aa_color_mapping(self):
        return {
            # Hydrophobic (green)
            "A": "green",
            "V": "green",
            "L": "green",
            "I": "green",
            "M": "green",
            "P": "cyan",
            # Aromatic
            "F": "magenta",
            "W": "magenta",
            "Y": "magenta",
            # Polar uncharged (blue)
            "N": "skyblue",
            "Q": "skyblue",
            "S": "skyblue",
            "T": "skyblue",
            # Negatively charged (red)
            "D": "red",
            "E": "red",
            # Positively charged (orange)
            "K": "blue",
            "R": "blue",
            "H": "orange",
            # Special cases
            "C": "gold",  # Often stands out (e.g., for disulfide bonds)
            "G": "gray",  # Smallest amino acid, neutral
            "X": "black",
        }


def _point_size(df, size, sizedef="area"):
    if sizedef == "area":
        def_size = 40
        def_max = 360
    else:
        def_size = 10
        def_max = 30
    # Define point size
    if size is None:
        size_values = def_size
    else:
        # Normalize size values (e.g., scale to [6, 30])
        min_size, max_size = def_size, def_max
        sizes = df[size].fillna(def_size)
        if sizes.max() == sizes.min():
            size_values = sizes
        else:
            norm_sizes = (sizes - sizes.min()) / (sizes.max() - sizes.min())
            size_values = norm_sizes * (max_size - min_size) + min_size
    return size_values


def _dual_marker_plotly(
    df,
    x="z-score",
    y="ali-length",
    color=None,
    marker=None,
    hovertext=None,
    title="Scatter Plot with Color and Marker Type",
    size=None,
):
    # Create figure
    fig = go.Figure()

    # Get unique color and marker categories
    colors = "gray"
    if color:
        colors = df[color].unique()
    symbols = None
    if marker:
        symbols = df[marker].unique()

    # Define marker symbols (extend as needed)
    marker_symbols = [
        "circle",
        "square",
        "diamond",
        "cross",
        "x",
        "triangle-up",
        "triangle-down",
    ]
    symbol_map = {
        s: marker_symbols[i % len(marker_symbols)] for i, s in enumerate(symbols)
    }

    # Define colors using Plotly default palette
    color_palette = [
        "#636EFA",
        "#00CC96",
        "#FFA15A",
        "#19D3F3",
        "#FF6692",
        "#B6E880",
        "#FF97FF",
        "#FECB52",
    ]
    color_map = {c: color_palette[i % len(color_palette)] for i, c in enumerate(colors)}

    # Create one trace per (color, symbol) combination
    for c in colors:
        for s in symbols:
            sub_df = df[(df[color] == c) & (df[marker] == s)]
            size_values = _point_size(sub_df, size, sizedef="diameter")
            if sub_df.empty:
                continue
            if hovertext is not None:
                fig.add_trace(
                    go.Scattergl(
                        x=sub_df[x],
                        y=sub_df[y],
                        mode="markers",
                        name=f"{c} | {s}",
                        marker=dict(
                            color=color_map[c],
                            symbol=symbol_map[s],
                            size=size_values,
                            line=dict(width=1, color="DarkSlateGrey"),
                        ),
                        text=sub_df[hovertext],
                        hoverinfo="text",
                    )
                )
            else:
                fig.add_trace(
                    go.Scattergl(
                        x=sub_df[x],
                        y=sub_df[y],
                        mode="markers",
                        name=f"{c} | {s}",
                        marker=dict(
                            color=color_map[c],
                            symbol=symbol_map[s],
                            size=size_values,
                            line=dict(width=1, color="DarkSlateGrey"),
                        ),
                    )
                )

    # Layout options
    fig.update_layout(
        title=title,
        xaxis_title=x,
        yaxis_title=y,
        width=900,
        height=600,
        template="plotly_white",
        legend=dict(
            x=1,
            y=0,
            # xanchor='right',
            # yanchor='bottom',
            bgcolor="rgba(255,255,255,0.6)",
            bordercolor="gray",
            borderwidth=1,
        ),
    )

    fig.show()


def read_TRUE(TRUEfile):
    "read hmmsearch results of target protein sequences against Pfam.35"
    TRUE_df = pd.read_csv(TRUEfile, sep="\t")
    TRUE_df.columns = TRUE_df.columns.str.lower()
    return TRUE_df


def diff_df(df1, df2, CLAN, TRUE_df, title="two-set comparison", plotfun="scatter"):
    """
    Compare two DataFrames by 'sbjct', summarize differences based on 'z-score', and visualize results.

    Parameters
    ----------
    df1 : pandas.DataFrame
        First dataset containing columns including 'sbjct', 'z-score', and 'ali-length'.

    df2 : pandas.DataFrame
        Second dataset with similar structure to df1.

    CLAN : object or None
        Optional clan information used for additional annotation in the plot. If None, no clan-based annotation is applied.

    TRUE_df : pandas.DataFrame
        DataFrame used for adding 'TRUE' annotations relative to the clan data.

    title : str, optional, default='two-set comparison'
        Title for the resulting plot.

    plotfun : str, optional, default='scatter'
        Type of plot to generate; accepts 'scatter' or other supported plotting methods (e.g., 'iscatter').

    Returns
    -------
    DALIFrame
        A custom DataFrame-like object containing merged data with added status labels and visualized plot.

    Details
    -------
    The function:
    - Extracts the highest 'z-score' row per 'sbjct' from each input DataFrame.
    - Performs an outer merge on 'sbjct' to align both datasets.
    - Classifies each 'sbjct' as:
        * 'both' if present in both datasets,
        * 'first-only' if only in df1,
        * 'second-only' if only in df2.
    - For overlapping entries, retains the row with the higher 'z-score' and corresponding 'ali-length'.
    - Optionally annotates with clan information via `TRUE_df`.
    - Generates a scatter plot or interactive scatter plot grouped by status and clan membership.

    Note
    ----
    - Requires `add_TRUE` and `DALIFrame` to be defined elsewhere in the codebase.
    - The function prints counts of each status category after annotation.
    """

    def get_status(row):
        if pd.notna(row["z-score_1"]) and pd.notna(row["z-score_2"]):
            return "both"
        elif pd.notna(row["z-score_1"]):
            return "first-only"
        else:
            return "second-only"

    def pick_values(row):
        if row["status"] == "second-only":
            return pd.Series([row["z-score_2"], row["ali-length_2"]])
        elif row["status"] == "first-only":
            return pd.Series([row["z-score_1"], row["ali-length_1"]])
        else:  # both
            if row["z-score_1"] >= row["z-score_2"]:
                return pd.Series([row["z-score_1"], row["ali-length_1"]])
            else:
                return pd.Series([row["z-score_2"], row["ali-length_2"]])

    def evaluate_two(df1, df2):
        # Step 1: Keep only the row with the highest z-score per sbjct in each df
        df1_max = df1.sort_values("z-score", ascending=False).drop_duplicates("sbjct")
        df2_max = df2.sort_values("z-score", ascending=False).drop_duplicates("sbjct")

        # Step 2: Outer join
        merged = pd.merge(
            df1_max, df2_max, on="sbjct", how="outer", suffixes=("_1", "_2")
        )

        # Step 3: Determine status
        merged["status"] = merged.apply(get_status, axis=1)

        # Step 4: Keep highest z-score and corresponding ali-length
        merged[["z-score", "ali-length"]] = merged.apply(pick_values, axis=1)

        return merged

    result = evaluate_two(df1, df2).drop_duplicates()[
        ["sbjct", "z-score", "ali-length", "status"]
    ]
    if CLAN is not None:
        # result is outer join so need to assign 'TRUE' (w.r.t. CLAN) if using it as marker in plot
        clean_df = add_TRUE(result, CLAN, TRUE_df).drop_duplicates()
        result = DALIFrame.from_df(clean_df)
        print(result["status"].value_counts())
        if plotfun == "scatter":
            result.scatter(color="status", marker="TRUE", title=title)
        else:
            result.iscatter(color="status", marker="TRUE", title=title)
    else:
        result = DALIFrame.from_df(result)
        if plotfun == "scatter":
            result.scatter(color="status", title=title)
        else:
            result.iscatter(color="status", title=title)
    return result


# plot structural dendrogram
# Use as: order=plot_ordered_occupancy_heatmap(df)
# pdist crashes if occupancy_df is large
#
def pileup_to_binary(pileup_str):
    return [0 if c == "." else 1 for c in pileup_str]


def plot_ordered_occupancy_heatmap(df, method="complete"):
    occupancy_df = df["dssp-pileup"].apply(pileup_to_binary)
    occupancy_array = pd.DataFrame(occupancy_df.tolist()).fillna(0).astype(int)

    # Compute linkage matrix using average linkage and Hamming distance (suitable for binary)
    linkage_matrix = linkage(pdist(occupancy_array, metric="hamming"), method=method)
    # Create clustermap
    g = sns.clustermap(
        occupancy_array,
        row_linkage=linkage_matrix,
        col_cluster=False,  # do not cluster columns (positions)
        cmap="Greys",
        figsize=(12, 8),
        # cbar=False,
        cbar_kws={
            "ticks": [],  # no ticks
            "label": "",  # no label
            "drawedges": False,  # optional: removes edges between color gradations
        },
    )
    plt.suptitle(f"Occupancy Heatmap with Row Clustering ({method} linkage)", y=1.02)
    plt.show()
    plt.close()
    return g.dendrogram_row.reordered_ind
