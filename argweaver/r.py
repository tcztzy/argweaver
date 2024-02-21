"""R interface for argweaver."""
import typing

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

from argweaver.utils import parse_region

if typing.TYPE_CHECKING:
    from os import PathLike
    from pathlib import Path
    from typing import Any, Dict, List, Literal, Optional, Tuple, Union

    import pandas

__all__ = ["plotTreesFromBed"]

_argweaver = importr("argweaver")


def plotTreesFromBed(
    bedfile: "Union[str, Path, PathLike]",
    *,
    i: "Union[int, Literal['max']]" = "max",
    s: "Optional[str]" = None,
) -> "Tuple[List[pandas.DataFrame], pandas.DataFrame]":
    """Plot trees from a bed file.

    See [](../../rapidocs/argweaver/plotTreesFromBed.md).

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
        chrom, start, end = parse_region(s)
        kwargs["chrom"] = chrom
        kwargs["start"] = start
        kwargs["end"] = end
    rv = _argweaver.plotTreesFromBed(str(bedfile), **kwargs)
    dfs = []
    for p in rv[0]:
        with (ro.default_converter + pandas2ri.converter).context():
            dfs.append(ro.conversion.get_conversion().rpy2py(p))
    with (ro.default_converter + pandas2ri.converter).context():
        x = ro.conversion.get_conversion().rpy2py(rv[1])
    return dfs, x
