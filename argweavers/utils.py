import io
import json
import re
import typing

import networkx as nx  # type: ignore
from Bio import Phylo

if typing.TYPE_CHECKING:
    from typing import Optional, Tuple, Union

    from Bio.Phylo.BaseTree import Clade, Tree

__all__ = ["parse_region", "parse_newick", "parse_nhx"]


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


def parse_nhx(nhx: str) -> dict:
    """Parse New Hampshire eXtended (NHX) format into a dictionary.

    Examples
    --------
    >>> parse_nhx("&&NHX:age=0.1")
    {'age': 0.1}
    """
    if str(nhx).startswith("&&NHX"):
        _, *kv = nhx.split(":")
        return json.loads(
            "{" + ",".join(['"' + pair.replace("=", '":') for pair in kv]) + "}"
        )
    else:
        return {}


def to_digraph(
    tree: "Union[Clade, Tree]", graph: "Optional[nx.DiGraph]" = None
) -> nx.DiGraph:
    """Convert tree to nx.DiGraph"""
    G = graph or nx.DiGraph()
    root = tree.root

    def add_edge(n1, n2, **kwargs):
        G.add_edge(n1, n2, **kwargs)

    G.add_node(
        root.name,
        pos=0,
        event="gene" if root.is_terminal() else "coal",
        **parse_nhx(root.comment),
    )
    if graph is None:
        G.nodes[root.name]["root"] = True
    for clade in root.clades:
        add_edge(clade.name, root.name)
        to_digraph(clade, G)
    return G


def parse_newick(newick: str) -> nx.DiGraph:
    """Parse a Newick format tree to a NetworkX DiGraph.

    Use this instead of Bio.Phylo.to_networkx for simplicity.
    Edge direction is from child to parent
    """
    tree = Phylo.read(io.StringIO(newick), "newick")
    for clade in tree.find_clades():
        clade.comment += ":local=true"
        if clade.is_terminal():
            clade.name = int(clade.name)
        else:
            clade.name = clade.confidence
            clade.confidence = None
    return to_digraph(tree)
