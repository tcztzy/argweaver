import importlib.util

import pytest

if importlib.util.find_spec("rpy2") is None:
    pytestmark = pytest.mark.skip(reason="rpy2 not installed")
else:
    from argweavers.r import plotTreesFromBed


def test_plot_trees(bedfile):
    with pytest.raises(ValueError):
        plotTreesFromBed(bedfile, s="not valid")
    rv = plotTreesFromBed(bedfile, s="chr:1000-2000")
    assert len(rv[0]) == len(rv[1])
