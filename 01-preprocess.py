#! /usr/bin/env python
import argparse
import glob
import json
import os
from collections import Counter, defaultdict

import pandas as pd


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
            counter = Counter([nt["molecule"] for nt in data["nucleotides"]])
            structure = data["quadruplexDotBracket"]["structure"].replace("-", "")
            chi = data["quadruplexDotBracket"]["chi"].replace("-", "")
            loop = data["quadruplexDotBracket"]["loop"].replace("-", "")

            for helix in data["helices"]:
                for quadruplex in helix["quadruplexes"]:
                    tetrads = quadruplex["tetrads"]
                    if len(tetrads) == 0:
                        continue

                    if any(
                        data["nucleotides"][nts[nt]]["shortName"] != "G"
                        for tetrad in tetrads
                        for nt in [
                            tetrad["nt1"],
                            tetrad["nt2"],
                            tetrad["nt3"],
                            tetrad["nt4"],
                        ]
                    ):
                        continue

                    qrs = [0] * len(nts)

                    for i, tetrad in enumerate(tetrads, 1):
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
                    assert len(qrs) == len(tetrads) * 4

                    key = (
                        counter.most_common(1)[0][0],
                        tuple(qrs),
                        tuple(description),
                        structure,
                        chi,
                        loop,
                    )
                    result[key].append(os.path.basename(json_path))

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
