from argweaver.plot import plot_tree


def test_plot_tree(bedfile):
    plot_tree(str(bedfile))
    assert True
