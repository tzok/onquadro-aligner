#! /usr/bin/env python
import argparse
import glob
import json
import os
from collections import Counter, defaultdict

import pandas as pd


def find_paths_of_length_n(graph, n):
    """
    Find all paths of length n in a directed acyclic graph.

    Args:
        graph (dict): A dictionary representing the graph as an adjacency list.
        n (int): The desired path length (number of edges).

    Returns:
        list: A list of all paths of length n.
    """
    all_paths = []

    def dfs(node, path, edges_traversed):
        # When we've traversed n-1 edges, we have a path of the desired length
        if edges_traversed == n - 1:
            all_paths.append(path.copy())
            return

        # Continue DFS if we haven't reached the desired length yet
        for neighbor in graph.get(node, []):
            path.append(neighbor)
            dfs(neighbor, path, edges_traversed + 1)
            path.pop()  # Backtrack

    # Start DFS from each node in the graph
    for start_node in graph:
        dfs(start_node, [start_node], 0)

    return all_paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "directory",
        help="Directory containing JSON files from ElTetrado",
    )
    parser.add_argument(
        "--output", help="Output file for results", default="onquadro-aligner.json"
    )
    args = parser.parse_args()
    result = defaultdict(list)

    for json_path in glob.iglob(f"{args.directory}/*.json"):
        with open(json_path) as f:
            data = json.load(f)

            if len(data["helices"]) == 0:
                continue

            nts = {nt["fullName"]: i for i, nt in enumerate(data["nucleotides"])}
            tetrads = {
                tetrad["id"]: tetrad
                for helix in data["helices"]
                for quadruplex in helix["quadruplexes"]
                for tetrad in quadruplex["tetrads"]
            }
            counter = Counter([nt["molecule"] for nt in data["nucleotides"]])

            pair_graph = defaultdict(list)
            for tetrad1, tetrad2 in [
                (pair["tetrad1"], pair["tetrad2"])
                for helix in data["helices"]
                for pair in helix["tetradPairs"]
            ]:
                pair_graph[tetrad1].append(tetrad2)

            n = 2
            while paths := find_paths_of_length_n(pair_graph, n):
                for path in paths:
                    flag = True
                    for tetrad_id in path:
                        tetrad = tetrads[tetrad_id]
                        if any(
                            data["nucleotides"][nts[nt]]["shortName"] != "G"
                            for nt in [
                                tetrad["nt1"],
                                tetrad["nt2"],
                                tetrad["nt3"],
                                tetrad["nt4"],
                            ]
                        ):
                            flag = False
                            break

                    if not flag:
                        continue

                    qrs = [0] * len(nts)

                    for i, tetrad_id in enumerate(path, 1):
                        tetrad = tetrads[tetrad_id]
                        qrs[nts[tetrad["nt1"]]] = i
                        qrs[nts[tetrad["nt2"]]] = i
                        qrs[nts[tetrad["nt3"]]] = i
                        qrs[nts[tetrad["nt4"]]] = i

                    current = []
                    description = []
                    for i, _ in enumerate(qrs):
                        if qrs[i] == 0:
                            current.append(data["nucleotides"][i]["shortName"])
                        else:
                            description.append(current.copy())
                            current = []
                    description.append(current.copy())
                    description = ["".join(d) for d in description]

                    qrs = [c for c in qrs if c != 0]
                    assert len(qrs) == len(path) * 4

                    structure = data["quadruplexDotBracket"]["structure"].replace(
                        "-", ""
                    )
                    chi = data["quadruplexDotBracket"]["chi"].replace("-", "")
                    loop = data["quadruplexDotBracket"]["loop"].replace("-", "")

                    key = (
                        counter.most_common(1)[0][0],
                        tuple(qrs),
                        tuple(description),
                        structure,
                        chi,
                        loop,
                    )
                    result[key].append(os.path.basename(json_path))
                n += 1

    rows = []
    for (molecule, qrs, description, structure, chi, loop), filenames in result.items():
        rows.append(
            {
                "molecule": molecule,
                "qrs": qrs,
                "description": description,
                "structure": structure,
                "chi": chi,
                "loop": loop,
                "filenames": filenames,
            }
        )

    df = pd.DataFrame(rows)
    df.to_json(args.output, orient="records")


if __name__ == "__main__":
    main()
