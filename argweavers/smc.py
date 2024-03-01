import gzip
import typing
from dataclasses import dataclass
from pathlib import Path

from argweavers.arg import ARG
from argweavers.utils import parse_newick

if typing.TYPE_CHECKING:
    from typing import List

__all__ = ["SMC"]


class SMC:
    """Sequential Markov Chain."""

    names: "List[str]"
    chrom: str
    start: int
    end: int
    trees: "List[Tree]"
    sprs: "List[SPR]"
    """Subtree Prune and Regrafts"""

    @dataclass
    class Tree:
        start: int
        end: int
        newick: str

        def __post_init__(self) -> None:
            self.G = parse_newick(self.newick)

    @dataclass
    class SPR:
        """Subtree Prune Re-graft."""

        pos: int
        """Position"""
        rnode: int
        """Recombination node"""
        rtime: float
        """Recombination time"""
        cnode: int
        """Coalescent node"""
        ctime: float
        """Coalescent time"""

    def __init__(self):
        self.trees = []
        self.sprs = []

    @classmethod
    def from_path(cls, path):
        path = Path(path)
        _open = gzip.open if path.suffix == ".gz" else open
        self = cls()
        with _open(path, "rt") as f:
            names_line = f.readline()
            tag, *self.names = names_line.split("\t")
            assert tag == "NAMES"
            region_line = f.readline()
            tag, self.chrom, start, end = region_line.split("\t")
            assert tag == "REGION"
            self.start, self.end = int(start), int(end)
            need = "TREE"
            for line in f:
                tag, *rest = line.split("\t")
                if tag == need == "TREE":
                    tree_start, tree_end, tree = rest
                    self.trees.append(cls.Tree(int(tree_start), int(tree_end), tree))
                    need = "SPR"
                elif tag == need == "SPR":
                    pos, recomb_node, recomb_time, coal_node, coal_time = rest
                    self.sprs.append(
                        (
                            int(pos),
                            int(recomb_node),
                            float(recomb_time),
                            int(coal_node),
                            float(coal_time),
                        )
                    )
                    need = "TREE"
                else:
                    raise ValueError(f"Unknown tag: {tag}")
            assert len(self.trees) == len(self.sprs) + 1
        return self

    def to_arg(self) -> ARG:
        raise NotImplementedError("WIP")
