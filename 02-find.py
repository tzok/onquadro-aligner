#! /usr/bin/env python
import argparse
import functools
import itertools
import json
import math
import multiprocessing
import os

import pandas as pd
from Bio import Align
from tqdm.contrib.concurrent import process_map


RESULT_COLUMNS = [
    "Files",
    "Tetrad count",
    "Molecule",
    "Tract distance",
    "Linker score",
    "Viability",
    "Loop lengths",
    "Topology",
    "Input sequence",
    "QRS",
]


VIABILITY_RANK = {
    "viable": 0,
    "marginal": 1,
    "unknown": 2,
    "n/a": 2,
    "not_viable": 3,
}


def _clamp(length):
    return min(length, 4)


def viability(l1, l2, l3, topology, n_tetrads):
    """Score whether a query fold is viable per two-tetrad G4 loop-length rules.

    Rules are derived from a manuscript studying two-tetrad DNA G4s with
    thymidine-only loops. The four observed topologies are applied to any
    tetrad count (lateral/diagonal loop energetics are largely tetrad-count
    independent). Propeller-containing topologies not among the four are
    treated as viable for >=3 tetrads (propeller loops are abundant in
    three-tetrad G4s) and not_viable for two tetrads (rare/unstable, except
    d+pd which is handled explicitly). Returns one of
    ``viable``/``marginal``/``not_viable``/``unknown``.
    """
    if n_tetrads < 2:
        return "n/a"

    a, b, c = _clamp(l1), _clamp(l2), _clamp(l3)

    if topology == "d+pd":
        return "viable" if a == 4 and c == 4 else "not_viable"

    if topology == "+l+l+l":
        if b == 3 and a in (2, 3, 4) and c in (2, 3, 4):
            if a == 4 and c == 4:
                return "not_viable"
            return "viable"
        if b == 3 and ((a == 1) ^ (c == 1)) and (
            a in (1, 2, 3, 4) and c in (1, 2, 3, 4)
        ):
            return "marginal"
        return "not_viable"

    if topology in ("-ld+l", "+ld-l"):
        if b == 4 and a in (3, 4) and c in (2, 3, 4):
            if a == 4 and c == 4:
                return "marginal"
            return "viable"
        if b == 4 and a == 2 and c in (2, 3, 4):
            return "marginal"
        return "not_viable"

    if topology == "-l-l-l":
        if b == 2 and a in (3, 4) and c in (3, 4):
            if a == 4 and c == 4:
                return "not_viable"
            return "marginal"
        return "not_viable"

    has_propeller = "p" in topology
    if has_propeller and topology not in ("d+pd",):
        if n_tetrads >= 3:
            if 1 in (a, b, c):
                return "not_viable"
            if a in (2, 3, 4) and b in (2, 3, 4) and c in (2, 3, 4):
                return "viable"
            return "marginal"
        return "not_viable"

    return "unknown"


def compute_query_loops(combination, gaps):
    """Return (l1, l2, l3) query loop lengths from the matched combination."""
    lengths = []
    for i, gap in enumerate(gaps):
        if gap + 1 >= len(combination):
            return None
        lengths.append(combination[gap + 1] - combination[gap] - 1)
    return tuple(lengths)


def assess_viability(obj, combination, n_tetrads):
    """Return (viability_label, loop_lengths_str, topology_str)."""
    loops = obj.get("loops")
    if not loops:
        return "n/a", "", ""
    gaps = loops["gaps"]
    topology = loops["topology"]
    lengths = compute_query_loops(combination, gaps)
    if lengths is None:
        return "not_viable", "", topology
    lengths_str = "-".join(str(length) for length in lengths)
    if any(length < 1 for length in lengths):
        return "not_viable", lengths_str, topology
    label = viability(lengths[0], lengths[1], lengths[2], topology, n_tetrads)
    return label, lengths_str, topology


@functools.cache
def combinations(iterable, r):
    return list(itertools.combinations(iterable, r))


@functools.cache
def describe(sequence, quadruplex_qrs, combination):
    qrs = [0] * len(sequence)

    for code, i in zip(quadruplex_qrs, combination):
        qrs[i] = code

    current = []
    description = []
    for i, _ in enumerate(qrs):
        if qrs[i] == 0:
            current.append(sequence[i])
        else:
            description.append(current.copy())
            current = []
    description.append(current.copy())
    return tuple(["".join(d) for d in description])


@functools.cache
def align(q, c):
    aligner = Align.PairwiseAligner(mode="global", scoring="blastn")
    # aligner.mode = "global"
    # aligner.open_gap_score = -1
    # aligner.extend_gap_score = -1
    return aligner.align(
        q.upper().replace("U", "T"),
        c.upper().replace("U", "T"),
    )


@functools.cache
def score_description(quadruplex, candidate):
    assert len(quadruplex) == len(candidate)

    score_tract = 0
    for i in range(1, len(quadruplex) - 1):
        if len(quadruplex[i]) == 0:
            score_tract += len(candidate[i])

    score_linkers = 0
    alignments = []
    for q, c in zip(quadruplex, candidate):
        if q and c:
            alignment = align(q, c)[0]
            score_linkers += int(alignment.score)
            _, a1, _, a2, _ = alignment.format("fasta").split("\n")
            alignments.append((a1, a2))
        elif q and not c:
            alignments.append((q, "-" * len(q)))
        elif not q and c:
            alignments.append(("-" * len(c), c))
        else:
            alignments.append(("", ""))

    return score_tract, score_linkers, alignments


def qrs_line(sequence_length, combination, qrs_chars):
    line = ["."] * sequence_length
    for i, char in zip(combination, qrs_chars):
        line[i] = char
    return "".join(line)


def serialize_g4composer_entry(
    name, sequence, structure, chi, sugar, orient, rise, twist, path, score_line=None
):
    fields = [
        ("name", name),
        ("sequence", sequence),
        ("structure", structure),
        ("chi", chi),
        ("sugar", sugar),
        ("orient", orient),
        ("rise", rise),
        ("twist", twist),
        ("path", path),
    ]
    body = "\n".join(f"{key:<11} {value}" for key, value in fields) + "\n"
    if score_line:
        return f"{score_line}\n{body}"
    return body


def write_g4composer_outputs(result, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    selected = []
    written = set()
    for row in result:
        g4c = row.get("g4c")
        if not g4c:
            continue
        base = row["template"]
        if base in written:
            continue
        written.add(base)
        selected.append(row)

    if not selected:
        return

    width = max(len(str(len(selected) - 1)), 1)

    for index, row in enumerate(selected):
        base = row["template"]
        g4c = row["g4c"]
        padded_score = str(index).zfill(width)

        sequence = row["Input sequence"]
        if row["Molecule"] == "RNA":
            seq = sequence.upper()
        else:
            seq = sequence.lower()

        structure = ["."] * len(sequence)
        for i in row["_combination"]:
            structure[i] = "^"
        structure = "".join(structure)

        score_line = (
            f"# tract_distance={row['Tract distance']} "
            f"linker_score={row['Linker score']} "
            f"viability={row['Viability']}"
        )
        if row["Loop lengths"]:
            score_line += f" loop_lengths={row['Loop lengths']}"
        if row["Topology"]:
            score_line += f" topology={row['Topology']}"
        score_line += f" template={base}"
        content = serialize_g4composer_entry(
            base,
            seq,
            structure,
            "." * len(sequence),
            "." * len(sequence),
            g4c["orient"],
            g4c["rise"],
            g4c["twist"],
            g4c["path"],
            score_line,
        )
        with open(os.path.join(output_dir, f"{padded_score}-{base}.inp"), "w") as f:
            f.write(content)


def process_single_item(args):
    """Process a single data item and return candidate results."""
    obj, sequence, g_indices, tetrad_count, list_all_quadruplex = args

    candidates = []
    files = " ".join(obj["filenames"]) if list_all_quadruplex else obj["filenames"][0]

    for n in tetrad_count:
        if len(obj["qrs"]) // 4 != n:
            continue

        key = (files, n)

        for combination in combinations(tuple(g_indices), n * 4):
            description = describe(sequence, tuple(obj["qrs"]), combination)
            score_tract, score_linkers, alignments = score_description(
                tuple(obj["description"]), description
            )
            candidates.append(
                (
                    key,
                    (
                        score_tract,
                        score_linkers,
                        obj,
                        combination,
                    ),
                )
            )

    return candidates


def match_quadruplexes(data, sequence, g_indices, tetrad_count, list_all_quadruplex):
    # Prepare arguments for parallel processing
    args_list = [
        (obj, sequence, g_indices, tetrad_count, list_all_quadruplex) for obj in data
    ]

    # Determine chunksize for better performance with large datasets
    chunksize = max(1, len(args_list) // (os.cpu_count() * 4))

    # Process in parallel with progress bar (uses all CPUs by default)
    all_candidates_lists = process_map(
        process_single_item,
        args_list,
        max_workers=os.cpu_count(),
        chunksize=chunksize,
        desc="Processing quadruplexes",
    )

    # Merge results
    best = {}
    for candidates in all_candidates_lists:
        for key, value in candidates:
            score_tract, score_linkers, obj, combination = value
            current = best.get(key, (math.inf, -math.inf, obj, []))

            if (score_tract, -score_linkers) < (current[0], -current[1]):
                best[key] = (score_tract, score_linkers, obj, combination)

    result = []

    for (files, n), (score_tract, score_linkers, obj, combination) in best.items():
        rendered_qrs = qrs_line(len(sequence), combination, tuple(obj["qrs_chars"]))
        v_label, v_lengths, v_topology = assess_viability(obj, combination, n)

        result.append(
            {
                "Files": files,
                "Tetrad count": n,
                "Molecule": obj["molecule"],
                "Tract distance": score_tract,
                "Linker score": score_linkers,
                "Viability": v_label,
                "Loop lengths": v_lengths,
                "Topology": v_topology,
                "Input sequence": sequence,
                "QRS": rendered_qrs,
                "g4c": obj.get("g4c"),
                "template": os.path.splitext(obj["filenames"][0])[0],
                "_combination": combination,
            }
        )

    return result


def print_results(result):
    if not result:
        print("No matches found.")
        return

    printed_qrs = set()
    print(result[0]["Input sequence"])
    for row in result:
        if row["QRS"] in printed_qrs:
            continue
        printed_qrs.add(row["QRS"])
        print(row["QRS"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        help="Input JSON file with preprocessed data",
        default="onquadro-aligner.json",
    )
    parser.add_argument(
        "--molecule",
        help="The quadruplex type to be analyzed against (default: any, choices: DNA, RNA, Other)",
    )
    parser.add_argument(
        "--tetrad-count",
        help="The number of tetrads to look for (default: maximum possible)",
        type=int,
    )
    parser.add_argument(
        "--list-all-quadruplex",
        help="Change the default behaviour of outputing only one non-redundant quadruplex ID",
        action="store_true",
    )
    parser.add_argument(
        "--output-csv",
        help="Write CSV results to this file",
    )
    parser.add_argument(
        "--g4composer-output-dir",
        help="Directory where g4composer .inp files are written (one per matched template)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=os.cpu_count(),
        help=f"Number of parallel workers (default: {os.cpu_count()} - all CPUs)",
    )
    parser.add_argument("sequence", help="The sequence to be analyzed")
    args = parser.parse_args()
    sequence = args.sequence.upper()

    with open(args.input) as f:
        data = json.load(f)

    if args.molecule:
        data = [d for d in data if d["molecule"] == args.molecule]

    g_indices = [i for i, c in enumerate(sequence) if c == "G"]
    max_tetrad_count = len(g_indices) // 4
    tetrad_count = (
        [args.tetrad_count]
        if args.tetrad_count
        else [max_tetrad_count]
    )

    result = match_quadruplexes(
        data,
        sequence,
        g_indices,
        tetrad_count,
        args.list_all_quadruplex,
    )
    result.sort(
        key=lambda row: (
            row["Tract distance"],
            -row["Linker score"],
            VIABILITY_RANK.get(row["Viability"], 2),
        )
    )
    print_results(result)

    if args.output_csv:
        pd.DataFrame(result, columns=RESULT_COLUMNS).to_csv(args.output_csv, index=False)

    if args.g4composer_output_dir:
        write_g4composer_outputs(result, args.g4composer_output_dir)


if __name__ == "__main__":
    main()
