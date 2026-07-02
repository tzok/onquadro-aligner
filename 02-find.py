#! /usr/bin/env python
import argparse
import functools
import itertools
import json
import math
import multiprocessing
import os
from dataclasses import dataclass
from typing import Optional, Union

import pandas as pd
from Bio import Align
from tqdm.contrib.concurrent import process_map


RESULT_COLUMNS = [
    "Files",
    "Tetrad count",
    "Molecule",
    "Tract distance",
    "Linker distance",
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

_PARALLEL = frozenset({"parallel"})
_CHAIR = frozenset({"chair"})
_BASKET = frozenset({"basket"})
_BASKET2 = frozenset({"basket2"})
_BASKET_BOTH = frozenset({"basket", "basket2"})
_HYBRID1 = frozenset({"hybrid1"})
_HYBRID3 = frozenset({"hybrid3"})
_HYBRID_BOTH = frozenset({"hybrid2", "hybrid3"})


def _clamp(length):
    return min(length, 4)


def parse_loop_types(topology):
    """Parse a topology string like ``"-ld+l"`` into ``["-l", "d", "+l"]``."""
    loops = []
    i = 0
    while i < len(topology):
        ch = topology[i]
        if ch == "d":
            loops.append("d")
            i += 1
        elif ch in "+-" and i + 1 < len(topology) and topology[i + 1] in "lp":
            loops.append(topology[i : i + 2])
            i += 2
        else:
            raise ValueError(f"Cannot parse topology at position {i}: {topology!r}")
    return loops


def loop_type_pattern(topology):
    """Return the (t1, t2, t3) type pattern, e.g. ``("p", "d", "l")``."""
    parsed = parse_loop_types(topology)
    return tuple("d" if t == "d" else t[1] for t in parsed)


def classify_topology(topology):
    """Map a loop progression to the set of topology classes it may belong to.

    Based on the canonical G4 topology classification (Silva et al.):
    8 groups, 26 topologies. Some patterns are ambiguous between two classes
    (e.g. (l,d,l) can be basket or basket2); in that case both are returned.
    Patterns not observed experimentally (e.g. (d,p,l)) return an empty set.
    """
    pattern = loop_type_pattern(topology)
    mapping = {
        ("p", "p", "p"): _PARALLEL,
        ("l", "l", "l"): _CHAIR,
        ("l", "d", "l"): _BASKET_BOTH,
        ("d", "p", "d"): _BASKET_BOTH,
        ("l", "p", "l"): _BASKET_BOTH,
        ("p", "l", "p"): _BASKET_BOTH,
        ("p", "d", "p"): _BASKET_BOTH,
        ("p", "l", "l"): _HYBRID1,
        ("p", "d", "l"): _HYBRID_BOTH,
        ("l", "d", "p"): _HYBRID_BOTH,
        ("p", "p", "l"): _HYBRID_BOTH,
        ("l", "l", "p"): _HYBRID_BOTH,
        ("l", "p", "p"): frozenset({"hybrid4"}),
    }
    return mapping.get(pattern, set())


@dataclass(frozen=True)
class _Rule:
    l1: Optional[int]
    l2: Optional[int]
    l3: Optional[int]
    cls: Optional[frozenset]
    n_tetrads: Optional[Union[int, tuple]]
    decision: str

    def matches(self, a, b, c, classes, n):
        if not _match(self.l1, a):
            return False
        if not _match(self.l2, b):
            return False
        if not _match(self.l3, c):
            return False
        if self.cls is not None:
            if not (self.cls & classes):
                return False
        else:
            if not classes:
                return False
        if not _match(self.n_tetrads, n):
            return False
        return True


def _match(condition, value):
    if condition is None:
        return True
    if isinstance(condition, tuple):
        op, val = condition
        if op == "ge":
            return value >= val
        if op == "ne":
            return value != val
        return False
    return value == condition


_RULES = [
    # --- Section A: Table 3 entries (l1 == l3) ---
    # (1, any, 1)
    _Rule(1, None, 1, _PARALLEL, ("ge", 3), "viable"),
    _Rule(1, None, 1, _PARALLEL, 2, "marginal"),
    _Rule(1, None, 1, None, None, "not_viable"),
    # (2, 1, 2)
    _Rule(2, 1, 2, _PARALLEL, ("ge", 3), "viable"),
    _Rule(2, 1, 2, _PARALLEL, 2, "marginal"),
    _Rule(2, 1, 2, None, None, "not_viable"),
    # (2, 2, 2)
    _Rule(2, 2, 2, _PARALLEL, ("ge", 3), "viable"),
    _Rule(2, 2, 2, _PARALLEL, 2, "marginal"),
    _Rule(2, 2, 2, _CHAIR, None, "marginal"),
    _Rule(2, 2, 2, None, None, "not_viable"),
    # (2, 3, 2)
    _Rule(2, 3, 2, _CHAIR, ("ne", 3), "viable"),
    _Rule(2, 3, 2, _CHAIR, 3, "marginal"),
    _Rule(2, 3, 2, None, None, "not_viable"),
    # (2, 4, 2)
    _Rule(2, 4, 2, _BASKET, ("ne", 3), "viable"),
    _Rule(2, 4, 2, _BASKET, 3, "marginal"),
    _Rule(2, 4, 2, None, None, "not_viable"),
    # (3, 1, 3)
    _Rule(3, 1, 3, _PARALLEL, ("ge", 3), "viable"),
    _Rule(3, 1, 3, _PARALLEL, 2, "marginal"),
    _Rule(3, 1, 3, None, None, "not_viable"),
    # (3, 2, 3)
    _Rule(3, 2, 3, _HYBRID3, ("ge", 3), "viable"),
    _Rule(3, 2, 3, _HYBRID3, 2, "marginal"),
    _Rule(3, 2, 3, None, None, "not_viable"),
    # (3, 3, 3)
    _Rule(3, 3, 3, _HYBRID1, ("ge", 3), "viable"),
    _Rule(3, 3, 3, _HYBRID1, 2, "marginal"),
    _Rule(3, 3, 3, _HYBRID3, None, "marginal"),
    _Rule(3, 3, 3, _CHAIR, None, "marginal"),
    _Rule(3, 3, 3, _BASKET2, None, "marginal"),
    _Rule(3, 3, 3, None, None, "not_viable"),
    # (3, 4, 3)
    _Rule(3, 4, 3, _BASKET, ("ne", 3), "viable"),
    _Rule(3, 4, 3, _BASKET, 3, "marginal"),
    _Rule(3, 4, 3, _BASKET2, None, "marginal"),
    _Rule(3, 4, 3, None, None, "not_viable"),
    # (4, 1, 4)
    _Rule(4, 1, 4, _BASKET, None, "marginal"),
    _Rule(4, 1, 4, None, None, "not_viable"),
    # (4, 4, 4)
    _Rule(4, 4, 4, _BASKET, None, "marginal"),
    _Rule(4, 4, 4, _CHAIR, None, "marginal"),
    _Rule(4, 4, 4, None, None, "not_viable"),
    # --- Section B: Two-tetrad rules (l1 != l3, n_tetrads=2) ---
    _Rule(4, 4, 4, _BASKET_BOTH, 2, "marginal"),
    _Rule(("ge", 3), 4, ("ge", 2), _BASKET_BOTH, 2, "viable"),
    _Rule(2, 4, ("ge", 2), _BASKET_BOTH, 2, "marginal"),
    _Rule(None, None, None, _BASKET_BOTH, 2, "not_viable"),
    _Rule(4, 2, 4, _CHAIR, 2, "not_viable"),
    _Rule(("ge", 3), 2, ("ge", 3), _CHAIR, 2, "marginal"),
    _Rule(None, None, None, _CHAIR, 2, "not_viable"),
    # --- Section C: Propeller fallback ---
    _Rule(("ge", 2), ("ge", 2), ("ge", 2), _PARALLEL, 2, "marginal"),
    _Rule(("ge", 2), ("ge", 2), ("ge", 2), _PARALLEL, ("ne", 2), "viable"),
    _Rule(None, None, None, _PARALLEL, None, "not_viable"),
    # --- Section D: Non-propeller fallback ---
    _Rule(("ge", 2), ("ge", 2), ("ge", 2), None, None, "marginal"),
]


def viability(l1, l2, l3, topology, n_tetrads):
    """Score whether a query fold is viable per gkag435 topological rules.

    Implements hard geometric constraints, then a single rules table
    (first-match-wins) covering Table 3 lookup, two-tetrad paper rules,
    tetrad-count adjustments, and general fallbacks. Returns one of
    ``viable``/``marginal``/``not_viable``/``unknown``/``n/a``.
    """
    if n_tetrads < 2:
        return "n/a"

    a, b, c = _clamp(l1), _clamp(l2), _clamp(l3)
    types = loop_type_pattern(topology)
    classes = classify_topology(topology)

    # --- Hard constraints (geometric impossibilities) ---
    for length, loop_type in zip((l1, l2, l3), types):
        if length <= 1 and loop_type != "p":
            return "not_viable"
    for length, loop_type in zip((l1, l2, l3), types):
        if loop_type == "d" and length < 4:
            return "not_viable"
    if n_tetrads >= 4 and any(length <= 1 for length in (l1, l2, l3)):
        return "not_viable"

    # --- Two-1-nt-loops rule: topology must be parallel ---
    short_count = sum(1 for length in (l1, l2, l3) if length <= 1)
    if short_count >= 2:
        if "parallel" in classes:
            return "viable"
        return "not_viable"

    # --- d+pd / d-pd explicit rule ---
    if topology in ("d+pd", "d-pd"):
        return "viable" if a == 4 and c == 4 else "not_viable"

    # --- Unrecognized topology (e.g. dpl — not observed) ---
    if not classes:
        return "unknown"

    # --- Rules table (first-match-wins) ---
    for rule in _RULES:
        if rule.matches(a, b, c, classes, n_tetrads):
            return rule.decision

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
    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = 0
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
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

    linker_distance = 0
    alignments = []
    for q, c in zip(quadruplex, candidate):
        if q and c:
            alignment = align(q, c)[0]
            linker_distance += -int(alignment.score)
            _, a1, _, a2, _ = alignment.format("fasta").split("\n")
            alignments.append((a1, a2))
        elif q and not c:
            linker_distance += len(q)
            alignments.append((q, "-" * len(q)))
        elif not q and c:
            linker_distance += len(c)
            alignments.append(("-" * len(c), c))
        else:
            alignments.append(("", ""))

    return score_tract, linker_distance, alignments


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
            f"linker_distance={row['Linker distance']} "
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
            score_tract, linker_distance, alignments = score_description(
                tuple(obj["description"]), description
            )
            candidates.append(
                (
                    key,
                    (
                        score_tract,
                        linker_distance,
                        obj,
                        combination,
                    ),
                )
            )

    return candidates


def match_quadruplexes(
    data, sequence, g_indices, tetrad_count, list_all_quadruplex, parallel=True
):
    # Prepare arguments for parallel processing
    args_list = [
        (obj, sequence, g_indices, tetrad_count, list_all_quadruplex) for obj in data
    ]

    if parallel:
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
    else:
        all_candidates_lists = [
            process_single_item(args) for args in args_list
        ]

    # Merge results — both metrics lower-is-better
    best = {}
    for candidates in all_candidates_lists:
        for key, value in candidates:
            score_tract, linker_distance, obj, combination = value
            current = best.get(key, (math.inf, math.inf, obj, []))

            if (score_tract, linker_distance) < (current[0], current[1]):
                best[key] = (score_tract, linker_distance, obj, combination)

    result = []

    for (files, n), (score_tract, linker_distance, obj, combination) in best.items():
        rendered_qrs = qrs_line(len(sequence), combination, tuple(obj["qrs_chars"]))
        v_label, v_lengths, v_topology = assess_viability(obj, combination, n)

        result.append(
            {
                "Files": files,
                "Tetrad count": n,
                "Molecule": obj["molecule"],
                "Tract distance": score_tract,
                "Linker distance": linker_distance,
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
            0
            if row["Tract distance"] == 0 and row["Linker distance"] == 0
            else 1,
            VIABILITY_RANK.get(row["Viability"], 2),
            row["Tract distance"],
            row["Linker distance"],
        )
    )
    print_results(result)

    if args.output_csv:
        pd.DataFrame(result, columns=RESULT_COLUMNS).to_csv(args.output_csv, index=False)

    if args.g4composer_output_dir:
        write_g4composer_outputs(result, args.g4composer_output_dir)


if __name__ == "__main__":
    main()
