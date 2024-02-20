import pytest

from argweaver.plot import plot_trees


def test_plot_trees(bedfile):
    with pytest.raises(ValueError):
        plot_trees(bedfile, s="not valid")
    rv = plot_trees(bedfile, s="chr:1000-2000")
    assert len(rv[0]) == len(rv[1])
