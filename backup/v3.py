#! /usr/bin/env python
import argparse
import functools
import glob
import itertools
import json
import math
import os
import re
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np


def read_quadruplexes(directory):
    result = []

    for path in glob.iglob(f"{directory}/*.json"):
        with open(path) as f:
            data = json.load(f)
            data["quadruplexDotBracket"]["filename"] = os.path.basename(path)
            if len(data["helices"]) > 0:
                data["quadruplexDotBracket"]["sequence"] = data["quadruplexDotBracket"][
                    "sequence"
                ].replace("-", "")
                data["quadruplexDotBracket"]["structure"] = data[
                    "quadruplexDotBracket"
                ]["structure"].replace("-", "")
                data["quadruplexDotBracket"]["chi"] = data["quadruplexDotBracket"][
                    "chi"
                ].replace("-", "")
                data["quadruplexDotBracket"]["loop"] = data["quadruplexDotBracket"][
                    "loop"
                ].replace("-", "")
                result.append(data["quadruplexDotBracket"])

    return result


def generate_tetrad_combinations(sequence):
    def is_valid_tetrad(x):
        return x[1] - x[0] > 1 and x[2] - x[1] > 1 and x[3] - x[2] > 1

    def is_valid_combination(x):
        used = set()
        for i in x:
            used.update(i)
        return len(used) == len(x) * 4

    indices = [i for i, x in enumerate(sequence) if x == "G"]
    tetrad_combinations = list(
        filter(is_valid_tetrad, itertools.combinations(indices, 4))
    )

    max_tetrads = len([x for x in sequence if x == "G"]) // 4

    return [
        x
        for i in range(max_tetrads, 1, -1)
        for x in itertools.combinations(tetrad_combinations, i)
        if is_valid_combination(x)
    ]


@functools.cache
def group_tetrad_combination(sequence, combination):
    s = list(" " * len(sequence))
    for i, tetrad in enumerate(combination):
        for index in tetrad:
            s[index] = chr(ord("q") + i)
    s = "".join(s)
    return s.strip().split(), [len(s) - len(s.lstrip())] + [
        len(x) for x in re.findall(r" +", s.strip())
    ] + [len(s) - len(s.rstrip())]


@functools.cache
def group_quadruplex(quadruplex):
    s = "".join(
        map(
            lambda c: c if "q" <= c <= "z" else " ",
            quadruplex.lower(),
        )
    )
    return s.strip().split(), [len(s) - len(s.lstrip())] + [
        len(x) for x in re.findall(r" +", s.strip())
    ] + [len(s) - len(s.rstrip())]


@functools.cache
def tetrad_combination_to_list(sequence, combination):
    # ((2, 5, 7, 11), (3, 6, 8, 19), (14, 20, 22, 24))
    # CUGGUGGGGCAGCAGCAAAGGGGAG
    # ..qr.qrqr..q..s....rs.s.s
    # [2, 0, 1, 0, 0, 0, 2, 2, 4, 0, 1, 1, 0]
    s = [False] * len(sequence)
    for tetrad in combination:
        for index in tetrad:
            s[index] = True
    count = 0
    result = []
    for b in s:
        if b:
            result.append(count)
            count = 0
        else:
            count += 1
    result.append(count)
    return result


@functools.cache
def quadruplex_to_list(structure):
    # GGGTGGGTTGGGTTGGG
    # QRS.qrs..QRS..qrs
    # [0, 0, 0, 1, 0, 0, 2, 0, 0, 2, 0, 0, 0]
    structure = structure.lower().replace("-", "")
    s = [False] * len(structure)
    for i, c in enumerate(structure):
        if "q" <= c <= "z":
            s[i] = True
    count = 0
    result = []
    for b in s:
        if b:
            result.append(count)
            count = 0
        else:
            count += 1
    result.append(count)
    return result


def similarity_score(list1, list2):
    score = 0
    for i in range(len(list1)):
        if list1[i] == 0 and list2[i] == 0:
            score += 2  # High reward for matching zeros
        elif list1[i] == 0 or list2[i] == 0:
            score -= 1  # Penalty for mismatched zeros
        else:
            # For non-zero values, reward similarity
            diff = abs(list1[i] - list2[i])
            max_val = max(list1[i], list2[i])
            score += 1 - (diff / max_val)

    return score / len(list1)  # Normalize


def zero_f1_score(list1, list2):
    # Count matching and mismatched zeros
    true_positives = sum(
        1 for i in range(len(list1)) if list1[i] == 0 and list2[i] == 0
    )
    false_positives = sum(
        1 for i in range(len(list1)) if list1[i] != 0 and list2[i] == 0
    )
    false_negatives = sum(
        1 for i in range(len(list1)) if list1[i] == 0 and list2[i] != 0
    )

    precision = (
        true_positives / (true_positives + false_positives)
        if (true_positives + false_positives) > 0
        else 0
    )
    recall = (
        true_positives / (true_positives + false_negatives)
        if (true_positives + false_negatives) > 0
        else 0
    )

    if precision + recall == 0:
        return 0

    return 2 * (precision * recall) / (precision + recall)


def weighted_distance(list1, list2):
    distance = 0

    for i in range(len(list1)):
        # Apply higher weight when zeros are involved
        weight = 3 if (list1[i] == 0 or list2[i] == 0) else 1
        distance += weight * (list1[i] - list2[i]) ** 2

    return 1 / (1 + math.sqrt(distance))  # Convert to similarity (0-1)


def score_combination(sequence, combination, quadruplex):
    l1 = tetrad_combination_to_list(sequence, combination)
    l2 = quadruplex_to_list(quadruplex["structure"])

    if len(l1) != len(l2):
        return False, (np.inf, 0, 0)

    g1 = group_tetrad_combination(sequence, combination)
    g2 = group_quadruplex(quadruplex["structure"])

    if "".join(g1[0]) != "".join(g2[0]):
        return False, (np.inf, 0, 0)

    return (
        True,
        (
            zero_f1_score(l1[1:-1], l2[1:-1]),
            weighted_distance(l1[1:-1], l2[1:-1]),
            similarity_score(l1[1:-1], l2[1:-1]),
        ),
    )

    group1 = group_tetrad_combination(sequence, combination)
    group2 = group_quadruplex(quadruplex["structure"])

    if "".join(group1[0]) == "".join(group2[0]):
        return True, np.linalg.norm(np.array(group1[1]) - np.array(group2[1]))

    return False, np.inf


def format_output(sequence, result):
    def process_group(group, seq):
        l1 = map(len, group[0])
        l2 = group[1]

        s1 = []
        begin = 0
        for i, j in itertools.zip_longest(l1, l2, fillvalue=0):
            s1.append(seq[begin : begin + j])
            s1.append(seq[begin + j : begin + j + i])
            begin += i + j
        return s1

    combination, quadruplex_group, _ = result
    group1 = group_tetrad_combination(sequence, combination)
    group2 = group_quadruplex(quadruplex_group[0]["structure"])

    s1 = []
    s2 = []
    s3 = []
    s4 = []
    s5 = []

    for g1, g2, g3, g4, g5 in zip(
        process_group(group1, sequence),
        process_group(group2, quadruplex_group[0]["sequence"]),
        process_group(group2, quadruplex_group[0]["structure"]),
        process_group(group2, quadruplex_group[0]["chi"]),
        process_group(group2, quadruplex_group[0]["loop"]),
    ):
        l = max(len(g1), len(g2))
        s1.append(f"{g1:->{l}}")
        s2.append(f"{g2:->{l}}")
        s3.append(f"{g3:->{l}}")
        s4.append(f"{g4:->{l}}")
        s5.append(f"{g5:->{l}}")

    result = f"Input:  {''.join(s1)}\n"
    result += f"Match:  {''.join(s2)}\n"
    result += f"        {''.join(s3)}\n"
    result += f"        {''.join(s4)}\n"
    result += f"        {''.join(s5)}\n"
    result += f"Source: {','.join([q['filename'] for q in quadruplex_group])}\n"
    return result


def process_combination(args):
    """Process a single combination against all quadruplex groups."""
    sequence, combination, unique_quadruplexes = args
    results = []

    for _, quadruplex_group in unique_quadruplexes.items():
        is_match, scores = score_combination(sequence, combination, quadruplex_group[0])
        if is_match:
            results.append((combination, quadruplex_group, scores, len(combination)))

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--directory",
        help="The directory containing the quadruplex JSON files",
        required=True,
    )
    parser.add_argument("sequence", help="The sequence to be analyzed")
    args = parser.parse_args()

    if not args.sequence:
        parser.print_help()
        sys.exit(1)

    tetrad_combinations = generate_tetrad_combinations(args.sequence)
    quadruplexes = read_quadruplexes(args.directory)

    unique_quadruplexes = defaultdict(list)
    for quadruplex in quadruplexes:
        unique_quadruplexes[(quadruplex["sequence"], quadruplex["structure"])].append(
            quadruplex
        )

    # Prepare arguments for parallel processing
    process_args = [
        (args.sequence, combination, unique_quadruplexes)
        for combination in tetrad_combinations
    ]

    # Use multiprocessing to distribute the workload
    results = defaultdict(list)
    with Pool() as pool:  # Use default number of processes (CPU count)
        for batch_results in pool.map(process_combination, process_args):
            for result in batch_results:
                if result:  # Only process non-empty results
                    combination, quadruplex_group, scores, num_tetrads = result
                    results[num_tetrads].append((combination, quadruplex_group, scores))

    for key in results:
        results[key].sort(key=lambda x: x[2], reverse=True)

    for key, value in results.items():
        print(f"Top match for {key} tetrads:")
        print(format_output(args.sequence, value[0]))
        print()


if __name__ == "__main__":
    main()
