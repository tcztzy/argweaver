"""Plotting functions for ARGweavers."""
import re
import typing

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

from argweaver.r import plotTreesFromBed

if typing.TYPE_CHECKING:
    from os import PathLike
    from pathlib import Path
    from typing import Any, Dict, Literal, Optional, Union

__all__ = ["plot_trees"]


def plot_trees(
    bedfile: "Union[str, Path, PathLike]",
    *,
    i: "Union[int, Literal['max']]" = "max",
    s: "Optional[str]" = None,
):
    """Plot trees from a bed file.

    Parameters
    ----------

    bedfile : {py:obj}`str`, {py:obj}`pathlib.Path` or {py:obj}`os.PathLike`
        Path to the bed file.

    i : {py:obj}`int`, optional
        The MCMC iteration to use. Default is -1, which means all intervals.

    s : {py:obj}`str`, optional
        Subset of the file, in format {chrom}:{start}-{end}
    """
    kwargs: "Dict[str, Any]" = {"iter": i}
    if s is not None:
        mo = re.match(r"(\w+):(\d+)-(\d+)", s)
        if mo is None:
            raise ValueError(f"Invalid subset string: {s}")
        kwargs["chrom"] = mo.group(1)
        kwargs["start"] = int(mo.group(2))
        kwargs["end"] = int(mo.group(3))
    rv = plotTreesFromBed(str(bedfile), **kwargs)
    dfs = []
    for p in rv[0]:
        with (ro.default_converter + pandas2ri.converter).context():
            dfs.append(ro.conversion.get_conversion().rpy2py(p))
    with (ro.default_converter + pandas2ri.converter).context():
        x = ro.conversion.get_conversion().rpy2py(rv[1])
    return dfs, x
