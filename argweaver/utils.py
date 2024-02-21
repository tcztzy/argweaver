import re
import typing

if typing.TYPE_CHECKING:
    from typing import Tuple

__all__ = ["parse_region"]


def parse_region(region: str) -> "Tuple[str, int, int]":
    """Parse a region string.

    Start and end are 1-based and inclusive. End must be greater than or equal to start.

    Parameters
    ----------
    region : str
        A string in the format {chrom}:{start}-{end}

    Returns
    -------
    Tuple[str, int, int]
        A tuple of (chrom, start, end)

    Raises
    ------
    ValueError
        If the region string is invalid.
    """
    mo = re.match(r"(\w+):(\d+)-(\d+)", region)
    if mo is None:
        raise ValueError(f"Invalid region string: {region}")
    chrom, start, end = mo.group(1), int(mo.group(2)), int(mo.group(3))
    if end < start:
        raise ValueError(f"Invalid region string: {region}")
    return chrom, start, end
