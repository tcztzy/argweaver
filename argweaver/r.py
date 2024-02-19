"""R interface for argweaver."""
from rpy2.robjects.packages import importr

_argweaver = importr("argweaver")

plotTreesFromBed = _argweaver.plotTreesFromBed
