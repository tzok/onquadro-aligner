#! /usr/bin/env python
import argparse
import functools
import itertools
import json
import math

import pandas as pd
from Bio import Align
from tqdm import tqdm


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
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    return aligner.align(
        q.upper(),
        c.upper(),
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


def match_quadruplexes(data, sequence, g_indices, tetrad_count, list_all_quadruplex):
    best = {}

    for obj in tqdm(data):
        files = (
            " ".join(obj["filenames"]) if list_all_quadruplex else obj["filenames"][0]
        )

        for n in tetrad_count:
            if len(obj["qrs"]) // 4 != n:
                continue

            key = (files, n)

            for combination in combinations(tuple(g_indices), n * 4):
                description = describe(sequence, tuple(obj["qrs"]), combination)
                score_tract, score_linkers, alignments = score_description(
                    tuple(obj["description"]), description
                )
                current = best.get(key, (math.inf, math.inf, obj, []))

                if (score_tract, score_linkers) < (current[0], current[1]) or \
                        (score_tract == current[0] and score_linkers > current[1]):
                    best[key] = (score_tract, score_linkers, obj, alignments)

    result = []

    for (files, n), (score_tract, score_linkers, obj, alignments) in best.items():
        aligned_input = "G".join([a[1] for a in alignments])
        aligned_quadruplex = "G".join([a[0] for a in alignments])
        structure = list(obj["structure"])
        chi = list(obj["chi"])
        loop = list(obj["loop"])

        for i, c in enumerate(aligned_quadruplex):
            if c == "-":
                structure.insert(i, "-")
                chi.insert(i, "-")
                loop.insert(i, "-")

        aligned_quadruplex += (
            "\n" + "".join(structure) + "\n" + "".join(chi) + "\n" + "".join(loop)
        )
        result.append(
            {
                "Files": files,
                "Tetrad count": n,
                "Molecule": obj["molecule"],
                "Tract distance": score_tract,
                "Linker score": score_linkers,
                "Aligned input": aligned_input,
                "Aligned quadruplex": aligned_quadruplex,
            }
        )

        assert sequence == result[-1]["Aligned input"].replace("-", ""), (
            sequence,
            result[-1]["Aligned input"].replace("-", ""),
        )

    return result


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
        help="The number of tetrads to look for (default: from 2 to maximum possible)",
        type=int,
    )
    parser.add_argument(
        "--list-all-quadruplex",
        help="Change the default behaviour of outputing only one non-redundant quadruplex ID",
        action="store_true",
    )
    parser.add_argument(
        "-o", "--output", help="Output file for results", default="results.csv"
    )
    parser.add_argument("sequence", help="The sequence to be analyzed")
    args = parser.parse_args()

    with open(args.input) as f:
        data = json.load(f)

    if args.molecule:
        data = [d for d in data if d["molecule"] == args.molecule]

    g_indices = [i for i, c in enumerate(args.sequence) if c == "G"]
    tetrad_count = (
        [args.tetrad_count]
        if args.tetrad_count
        else list(range(2, len(g_indices) // 4 + 1))
    )

    result = match_quadruplexes(
        data,
        args.sequence,
        g_indices,
        tetrad_count,
        args.list_all_quadruplex,
    )
    df = pd.DataFrame(result)
    df.sort_values(["Tract distance", "Linker score"], ascending=[True, False], inplace=True)
    df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
