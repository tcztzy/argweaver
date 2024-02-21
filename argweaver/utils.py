import re
import typing

if typing.TYPE_CHECKING:
    from typing import Tuple


def parse_region(region: str) -> "Tuple[str, int, int]":
    mo = re.match(r"(\w+):(\d+)-(\d+)", region)
    if mo is None:
        raise ValueError(f"Invalid region string: {region}")
    chrom, start, end = mo.group(1), int(mo.group(2)), int(mo.group(3))
    if end < start:
        raise ValueError(f"Invalid region string: {region}")
    return chrom, start, end
