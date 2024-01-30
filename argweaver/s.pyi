from pathlib import Path
from typing import Union

import polars as pl

class Sites:
    chrom: str
    start: int
    end: int
    data: pl.DataFrame

def read_sites(filename: Union[str, Path]) -> Sites: ...
