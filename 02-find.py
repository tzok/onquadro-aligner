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
    "Input sequence",
    "QRS",
]


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

        result.append(
            {
                "Files": files,
                "Tetrad count": n,
                "Molecule": obj["molecule"],
                "Tract distance": score_tract,
                "Linker score": score_linkers,
                "Input sequence": sequence,
                "QRS": rendered_qrs,
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
    result.sort(key=lambda row: (row["Tract distance"], -row["Linker score"]))
    print_results(result)

    if args.output_csv:
        pd.DataFrame(result, columns=RESULT_COLUMNS).to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    main()
