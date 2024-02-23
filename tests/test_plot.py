import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison

from argweavers.plot import plot_tree


@image_comparison(baseline_images=["test_plot_tree"], extensions=["png"])
def test_plot_tree():
    fig, ax = plt.subplots()
    plot_tree("((A:1,B:1):1,C:2);", name="Test Plot", axes=ax, do_show=False)
