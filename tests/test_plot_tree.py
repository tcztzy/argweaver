from argweaver.plot import plot_trees


def test_plot_tree(bedfile):
    plot_trees(bedfile)
    assert True
