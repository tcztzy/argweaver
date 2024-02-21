"""Plotting functions for ARGweavers."""
import typing
from io import StringIO

from Bio import Phylo

if typing.TYPE_CHECKING:
    pass

__all__ = ["plot_tree"]


def plot_tree(newick: str, *, name="", **kwargs):
    """Plot a tree

    Parameters
    ----------
    newick : {py:obj}`str`
        A newick string containing a tree.
    name : {py:obj}`str`, optional
        Name of the tree.
    **kwargs
        Additional keyword arguments to pass to {py:obj}`Bio.Phylo.draw`.
    """
    tree = Phylo.read(StringIO(newick), "newick")
    if name:
        tree.name = name
    tree.ladderize()
    Phylo.draw(tree, **kwargs)
