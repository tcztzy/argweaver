import gzip
from pathlib import Path

import pandas as pd

from argweaver.utils import parse_region

__all__ = ["read_bed"]


def read_bed(filename):
    """Read a bed file.

    Parameters
    ----------
    filename : {py:obj}`str`
        Path to the bed file.

    Returns
    -------
    {py:obj}`pandas.DataFrame`
    """
    path = Path(filename)
    _open = gzip.open if path.suffix == ".gz" else open
    return pd.read_table(
        _open(path, "rt"), header=None, names=["chrom", "start", "end", "iter", "tree"]
    )


def subset_bed(data, region):
    """Subset a bed file.

    Parameters
    ----------
    bedfile : {py:obj}`str`
        Path to the bed file.
    region : {py:obj}`str`
        Subset of the file, in format {chrom}:{start}-{end}

    Returns
    -------
    {py:obj}`pandas.DataFrame`
    """
    chrom, start, end = parse_region(region)
    intervals = pd.IntervalIndex.from_arrays(data["start"], data["end"], closed="left")
    interval = pd.Interval(start, end, closed="left")
    return data[(data["chrom"] == chrom) & intervals.overlaps(interval)]
