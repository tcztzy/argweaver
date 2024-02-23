from pathlib import Path
from typing import List, Optional, Union

import polars as pl

__all__: list = ["Sites", "read_sites", "smc2bed"]

class Sites:
    """A class to represent a set of sites in a region.

    Attributes
    ----------
    chrom : str
        The chromosome name.
    start : int
        The start position of the region.
    end : int
        The end position of the region.
    data : polars.DataFrame

    Notes
    -----
    Sites instance can be sliced to get a new Sites instance with the same
    chrom, but new start, and new end and a subset of the data.

    Examples
    --------
    >>> sites = read_sites("examples/sim1/sim1.sites")
    >>> assert sites.data.shape == (170, 9)
    >>> sites_copy = sites[405:]
    >>> assert sites_copy.data.shape == (169, 9)
    >>> assert sites_copy.start == 405
    """

    chrom: str
    start: int
    end: int
    data: pl.DataFrame

    def __getitem__(self, s: slice) -> Sites: ...

def read_sites(filename: Union[str, Path]) -> Sites: ...
def smc2bed(args: Optional[List[str]]) -> None: ...
