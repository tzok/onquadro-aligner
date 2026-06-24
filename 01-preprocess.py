#! /usr/bin/env python
import argparse
import glob
import json
import math
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


def remap_path_entry_to_clockwise(entry):
    label = "".join(ch for ch in entry if ch.isalpha())
    column = "".join(ch for ch in entry if ch.isdigit())
    if column == "1":
        return entry
    clockwise_column = {"2": "4", "3": "3", "4": "2"}.get(column)
    if clockwise_column is None:
        raise ValueError(f"Unsupported g4composer path column: {entry}")
    return f"{label}{clockwise_column}"


def format_number(value):
    if math.isnan(value):
        return "."
    rounded = round(value, 1)
    if math.isclose(rounded, round(rounded), abs_tol=1.0e-9):
        return str(int(round(rounded)))
    return f"{rounded:.1f}".rstrip("0").rstrip(".")


def tetrad_label(index):
    labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if index >= len(labels):
        raise ValueError(f"Cannot represent more than {len(labels)} tetrads")
    return labels[index]


def path_steps_between(start, end, adjacency, steps):
    queue = [(start, [], [start])]
    while queue:
        current, path_steps, visited = queue.pop(0)
        if current == end:
            return path_steps
        for neighbor in adjacency.get(current, []):
            if neighbor in visited:
                continue
            queue.append(
                (
                    neighbor,
                    path_steps + [steps[(current, neighbor)]],
                    visited + [neighbor],
                )
            )
    raise ValueError("Failed to derive g4composer build traversal between tetrads")


def compute_g4composer(quadruplex, tetrad_pairs, global_index):
    """Reconstruct g4composer orient/rise/twist/path from ElTetrado JSON data.

    Mirrors eltetrado.g4composer but operates purely on the JSON DTO so no 3D
    coordinates are required. The JSON tetrad-pair rise is already the unsigned
    magnitude that g4composer uses as signed-rise magnitude; interval signs come
    only from build-order traversal direction.
    """
    tetrads = quadruplex["tetrads"]
    if len(tetrads) < 2:
        return None

    polarities = quadruplex.get("tetradPolarities") or []
    if len(polarities) != len(tetrads) or any(p is None for p in polarities):
        return None

    tetrad_ids = [t["id"] for t in tetrads]
    tetrad_min_index = {
        t["id"]: min(
            global_index[t["nt1"]],
            global_index[t["nt2"]],
            global_index[t["nt3"]],
            global_index[t["nt4"]],
        )
        for t in tetrads
    }
    ordered_ids = sorted(tetrad_ids, key=lambda tid: tetrad_min_index[tid])
    id_to_label = {tid: tetrad_label(i) for i, tid in enumerate(ordered_ids)}
    label_to_id = {label: tid for tid, label in id_to_label.items()}

    path = quadruplex["path"]
    build_labels = []
    for entry in path:
        label = "".join(ch for ch in entry if ch.isalpha())
        if label not in build_labels:
            build_labels.append(label)
    build_order_ids = [label_to_id[label] for label in build_labels]

    polarity_by_id = {t["id"]: p for t, p in zip(tetrads, polarities)}
    orient = ";".join(
        f"{id_to_label[tid]}{'+' if polarity_by_id[tid] == 'clockwise' else '-'}"
        for tid in build_order_ids
    )

    steps = {}
    for tp in tetrad_pairs:
        a, b = tp["tetrad1"], tp["tetrad2"]
        rise = tp["rise"]
        twist = tp["twist"]
        steps[(a, b)] = (rise, twist)
        steps[(b, a)] = (-rise, -twist)
    adjacency = {}
    for a, b in steps:
        adjacency.setdefault(a, []).append(b)

    rise_values = []
    twist_values = []
    for start, end in zip(build_order_ids, build_order_ids[1:]):
        interval = path_steps_between(start, end, adjacency, steps)
        rise_values.append(format_number(sum(s[0] for s in interval)))
        twist_values.append(format_number(sum(s[1] for s in interval)))

    return {
        "orient": orient,
        "rise": ";".join(rise_values),
        "twist": ";".join(twist_values),
        "path": ";".join(remap_path_entry_to_clockwise(e) for e in path),
    }


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
        "--output", help="Output file for results", default="onquadro-aligner.json"
    )
    args = parser.parse_args()
    result = defaultdict(set)
    g4c_by_key = {}
    loops_by_key = {}

    for json_path in glob.iglob(f"{args.directory}/*.json"):
        with open(json_path) as f:
            data = json.load(f)

        if len(data["helices"]) == 0:
            continue

        nts = {nt["fullName"]: i for i, nt in enumerate(data["nucleotides"])}
        global_index = {nt["fullName"]: nt["index"] for nt in data["nucleotides"]}
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
                    tetrad_ids = {tetrad["id"] for tetrad in tetrads}
                    quadruplex_tetrad_pairs = [
                        tp
                        for tp in helix["tetradPairs"]
                        if tp["tetrad1"] in tetrad_ids
                        and tp["tetrad2"] in tetrad_ids
                    ]
                    g4c_by_key[key] = compute_g4composer(
                        quadruplex, quadruplex_tetrad_pairs, global_index
                    )
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
