#! /usr/bin/env python
import argparse
import glob
import json
import os
from collections import Counter, defaultdict

import pandas as pd
from tqdm import tqdm


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


def parse_g4composer_inp(path):
    """Read orient/rise/twist/path directly from a g4composer .inp file."""
    fields = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split(maxsplit=1)
            if len(parts) != 2:
                continue
            key, value = parts
            if key in ("orient", "rise", "twist", "path"):
                fields[key] = value
    if not {"orient", "rise", "twist", "path"}.issubset(fields):
        return None
    return fields


LOOP_TYPE_MAP = {
    "lateral+": "+l",
    "lateral-": "-l",
    "diagonal": "d",
    "propeller+": "+p",
    "propeller-": "-p",
}


def compute_loops_info(quadruplex, chain_positions, chain):
    """Reconstruct loop topology and gap indices from ElTetrado JSON data.

    Returns ``{"topology": str, "gaps": [k1, k2, k3]}`` or ``None``. The gap
    index ki is the position, in the chain-sorted list of all tetrad
    nucleotides, of the tetrad G immediately 5' of loop i. The query loop
    length is then ``combination[gaps[i]+1] - combination[gaps[i]] - 1``.
    """
    loops = quadruplex.get("loops") or []
    if len(loops) != 3 or any(loop.get("type") is None for loop in loops):
        return None

    try:
        topology = "".join(LOOP_TYPE_MAP[loop["type"]] for loop in loops)
    except KeyError:
        return None

    tetrad_nts = []
    for tetrad in quadruplex["tetrads"]:
        for k in range(1, 5):
            tetrad_nts.append(tetrad[f"nt{k}"])
    sorted_tetrad = sorted(
        tetrad_nts, key=lambda nt: chain_positions[chain][nt]
    )

    gaps = []
    for loop in loops:
        first_nt = loop["nucleotides"][0]
        loop_pos = chain_positions[chain][first_nt]
        gap_idx = None
        for j, tn in enumerate(sorted_tetrad):
            if chain_positions[chain][tn] < loop_pos:
                gap_idx = j
            else:
                break
        if gap_idx is None:
            return None
        gaps.append(gap_idx)

    return {"topology": topology, "gaps": gaps}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "directory",
        help="Directory containing JSON files from ElTetrado",
    )
    parser.add_argument(
        "--g4composer-dir",
        required=True,
        help="Directory containing g4composer .inp files paired (by base name) "
        "with the ElTetrado JSON files",
    )
    parser.add_argument(
        "--output", help="Output file for results", default="onquadro-aligner.json"
    )
    args = parser.parse_args()
    result = defaultdict(set)
    g4c_by_key = {}
    loops_by_key = {}

    inp_paths = sorted(glob.glob(f"{args.g4composer_dir}/*.inp"))
    for inp_path in tqdm(inp_paths, desc="Preprocessing"):
        base = os.path.splitext(os.path.basename(inp_path))[0]
        json_path = os.path.join(args.directory, f"{base}.json")
        if not os.path.exists(json_path):
            continue

        g4c = parse_g4composer_inp(inp_path)
        if g4c is None:
            continue

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
                if key not in g4c_by_key:
                    g4c_by_key[key] = g4c
                    loops_by_key[key] = compute_loops_info(
                        quadruplex, chain_positions, chain
                    )

    rows = []
    for (molecule, qrs, qrs_chars, description), filenames in result.items():
        rows.append(
            {
                "molecule": molecule,
                "qrs": qrs,
                "qrs_chars": qrs_chars,
                "description": description,
                "filenames": sorted(filenames),
                "g4c": g4c_by_key[(molecule, qrs, qrs_chars, description)],
                "loops": loops_by_key[(molecule, qrs, qrs_chars, description)],
            }
        )

    df = pd.DataFrame(rows)
    df.to_json(args.output, orient="records")


if __name__ == "__main__":
    main()
