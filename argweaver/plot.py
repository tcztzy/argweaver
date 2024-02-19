"""Plotting functions for ARGweavers."""
from argweaver.r import plotTreesFromBed


def plot_trees(bedfile: str):
    """Plot trees from a bed file.

    Parameters
    ----------

    bedfile : {py:obj}`str` or {py:obj}`pathlib.Path`
        Path to the bed file.
    """
    plotTreesFromBed(str(bedfile))
