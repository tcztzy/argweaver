"""R interface for argweaver."""
from rpy2.robjects.packages import importr

__all__ = ["plotTreesFromBed"]

_argweaver = importr("argweaver")

plotTreesFromBed = _argweaver.plotTreesFromBed
"""Plot trees from a bed file.

See [](../../rapidocs/argweaver/plotTreesFromBed.md).
"""
