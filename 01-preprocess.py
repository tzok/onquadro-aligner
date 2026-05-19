#! /usr/bin/env python
import argparse
import glob
import json
import os
from collections import Counter, defaultdict

import pandas as pd


def chain_id(full_name):
    return full_name.split(".")[0]


def describe(sequence, qrs):
    current = []
    description = []

    for nt, code in zip(sequence, qrs):
        if code == 0:
            current.append(nt)
        else:
            description.append("".join(current))
            current = []

    description.append("".join(current))
    return tuple(description)


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
    result = defaultdict(set)

    for json_path in glob.iglob(f"{args.directory}/*.json"):
        with open(json_path) as f:
            data = json.load(f)

        if len(data["helices"]) == 0:
            continue

        nts = {nt["fullName"]: i for i, nt in enumerate(data["nucleotides"])}
        structure = data["quadruplexDotBracket"]["structure"].replace("-", "")

        chains = defaultdict(list)
        for i, nucleotide in enumerate(data["nucleotides"]):
            chains[chain_id(nucleotide["fullName"])].append(i)

        chain_sequences = {}
        chain_positions = {}
        chain_molecules = {}
        for chain, indices in chains.items():
            chain_sequences[chain] = [data["nucleotides"][i]["shortName"] for i in indices]
            chain_positions[chain] = {
                data["nucleotides"][i]["fullName"]: j for j, i in enumerate(indices)
            }
            chain_molecules[chain] = Counter(
                data["nucleotides"][i]["molecule"] for i in indices
            ).most_common(1)[0][0]

        for helix in data["helices"]:
            for quadruplex in helix["quadruplexes"]:
                tetrads = quadruplex["tetrads"]
                if len(tetrads) == 0:
                    continue

                quadruplex_nts = [
                    nt
                    for tetrad in tetrads
                    for nt in [
                        tetrad["nt1"],
                        tetrad["nt2"],
                        tetrad["nt3"],
                        tetrad["nt4"],
                    ]
                ]
                quadruplex_chains = {chain_id(nt) for nt in quadruplex_nts}
                if len(quadruplex_chains) != 1:
                    continue

                chain = next(iter(quadruplex_chains))
                qrs = [0] * len(chain_sequences[chain])
                qrs_chars = ["."] * len(chain_sequences[chain])
                valid = True

                for i, tetrad in enumerate(tetrads, 1):
                    for nt in [
                        tetrad["nt1"],
                        tetrad["nt2"],
                        tetrad["nt3"],
                        tetrad["nt4"],
                    ]:
                        index = nts[nt]
                        if data["nucleotides"][index]["shortName"] != "G":
                            valid = False
                            break

                        local_index = chain_positions[chain][nt]
                        qrs[local_index] = i
                        qrs_chars[local_index] = structure[index]
                        assert qrs_chars[local_index].isalpha(), (
                            os.path.basename(json_path),
                            nt,
                            qrs_chars[local_index],
                        )

                    if not valid:
                        break

                if not valid:
                    continue

                description = describe(chain_sequences[chain], qrs)
                qrs = tuple(code for code in qrs if code != 0)
                qrs_chars = tuple(char for char in qrs_chars if char != ".")
                assert len(qrs) == len(tetrads) * 4
                assert len(qrs_chars) == len(tetrads) * 4

                key = (
                    chain_molecules[chain],
                    qrs,
                    qrs_chars,
                    description,
                )
                result[key].add(os.path.basename(json_path))

    rows = []
    for (molecule, qrs, qrs_chars, description), filenames in result.items():
        rows.append(
            {
                "molecule": molecule,
                "qrs": qrs,
                "qrs_chars": qrs_chars,
                "description": description,
                "filenames": sorted(filenames),
            }
        )

    df = pd.DataFrame(rows)
    df.to_json(args.output, orient="records")


if __name__ == "__main__":
    main()
