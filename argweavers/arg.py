import re
from pathlib import Path

import networkx as nx  # type: ignore
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import tskit  # type: ignore

__all__ = ["ARG"]


class ARG:
    G: nx.DiGraph

    @property
    def sequence_length(self) -> int:
        return self.G["sequence_length"]

    @classmethod
    def from_path(cls, path):
        """Construct an ARG from an ARGweaver file.

        The `.arg` file start with a header line with the start and end positions of the sequence.

        The rest of the file is a tab-separated table with the following columns:
        - name: the name of the node
        - event: the type of event (gene, recomb, coal)
        - age: the age of the node
        - pos: the position of the recombination event (if event is recomb)
        - parents: the parents of the node
        - children: the children of the node
        """
        path = Path(path)
        self = cls()
        with open(path) as f:
            line = next(f).strip()
            mo = re.match(r"^start=(\d+)\tend=(\d+)$", line)
            assert mo is not None

            # the "name" field can be a string. Force it to be so, in case it is just numbers
            df = pd.read_csv(
                f, header=0, sep="\t", dtype={"name": str, "parents": str}, index_col=0
            )
            for col in ("parents", "age"):
                if col not in df.columns:
                    raise ValueError(f"Column {col} not found in ARGweaver file")

            names = sorted(df.query("event == 'gene'").index)
            # Make an nx DiGraph so we can do a topological sort.
            G = nx.DiGraph()
            assert int(mo.group(1)) == 0
            G["sequence_length"] = int(mo.group(2))
            for child, row in df.iterrows():
                try:
                    child = int(child)
                except ValueError:
                    pass
                if row["event"] == "gene":
                    G.add_node(
                        names.index(child),
                        age=row["age"],
                        event=row["event"],
                        sample=child,
                    )
                elif row["event"] == "recomb":
                    G.add_node(
                        int(child), age=row["age"], event=row["event"], pos=row["pos"]
                    )
                else:
                    G.add_node(int(child), age=row["age"], event=row["event"])
                if isinstance(row["parents"], str):
                    for i, parent in enumerate(
                        [int(p) for p in row["parents"].split(",")]
                    ):
                        if child in names:
                            child = names.index(child)
                        G.add_edge(child, parent, parent=i)
                else:
                    G.add_node(child, root=True)
            self.G = G
        return self

    def to_ts(self) -> tskit.TreeSequence:
        tables = tskit.TableCollection(sequence_length=self.sequence_length)
        tables.nodes.metadata_schema = tskit.MetadataSchema.permissive_json()
        breakpoints = np.full(len(self.G), tables.sequence_length)
        aw_to_tsk_id = {}
        times = np.unique([a for _, a in self.G.nodes.data("age")])
        time_map = {time: 1 if time == 0 else 0 for time in times}
        min_time_diff = min(np.diff(times))
        epsilon = min_time_diff / 1e6
        try:
            for node in nx.lexicographical_topological_sort(self.G):
                record = self.G.nodes[node].copy()
                record["name"] = record.get("sample", str(node))
                age = record["age"]
                flags = 0
                # Sample nodes are marked as "gene" events
                if record["event"] == "gene":
                    flags = tskit.NODE_IS_SAMPLE
                    assert age == 0
                    time = age
                else:
                    time = age + time_map[age] * epsilon
                    # Argweaver allows age of parent and child to be the same, so we
                    # need to add epsilons to enforce parent_age > child_age
                    time_map[age] += 1
                tsk_id = tables.nodes.add_row(flags=flags, time=time, metadata=record)
                aw_to_tsk_id[node] = tsk_id
                if record["event"] == "recomb":
                    breakpoints[tsk_id] = record["pos"]
        except nx.exception.NetworkXUnfeasible:
            bad_edges = nx.find_cycle(self.G, orientation="original")
            raise nx.exception.NetworkXUnfeasible(
                f"Cycle found in ARGweaver graph: {bad_edges}"
            )

        L = tables.sequence_length
        for aw_node in self.G:
            child = aw_to_tsk_id[aw_node]
            parents = [
                aw_to_tsk_id[aw_parent]
                for _, aw_parent, _ in sorted(
                    self.G.edges(aw_node, data="parent"), key=lambda x: x[2]
                )
            ]
            if len(parents) == 1:
                tables.edges.add_row(0, L, parents[0], child)
            elif len(parents) == 2:
                # Recombination node.
                # Note that this uses the 1-RE-node convention
                x = breakpoints[child]
                tables.edges.add_row(0, x, parents[0], child)
                tables.edges.add_row(x, L, parents[1], child)
            else:
                assert len(parents) == 0
        tables.sort()
        ts = tables.tree_sequence()
        return ts.simplify(keep_unary=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    ARG.from_path(args.path).to_ts().dump(args.path.with_suffix(".trees"))
