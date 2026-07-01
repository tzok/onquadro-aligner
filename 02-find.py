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

    Based on Table 1 of gkag435. Some patterns are ambiguous between two
    classes (e.g. (l,l,p) can be hybrid3 or hybrid4); in that case both are
    returned.
    """
    pattern = loop_type_pattern(topology)
    mapping = {
        ("p", "p", "p"): {"parallel"},
        ("l", "l", "l"): {"chair"},
        ("l", "d", "l"): {"basket", "basket2"},
        ("d", "p", "d"): {"basket"},
        ("l", "p", "l"): {"basket2"},
        ("p", "d", "p"): {"hybrid1"},
        ("p", "l", "l"): {"hybrid1"},
        ("p", "d", "l"): {"hybrid2"},
        ("l", "d", "p"): {"hybrid2"},
        ("p", "p", "l"): {"hybrid2", "hybrid3"},
        ("l", "l", "p"): {"hybrid3", "hybrid4"},
        ("d", "p", "l"): {"hybrid4"},
        ("l", "p", "p"): {"hybrid4"},
    }
    return mapping.get(pattern, set())


# Table 3 from gkag435: (clamped L1, clamped L2, clamped L3) -> {class: percentage}
# "any" means any L2 value. L2=4 covers the 4-5 nt range (clamped).
TABLE3 = {
    (1, "any", 1): {"parallel": 100},
    (2, 1, 2): {"parallel": 100},
    (2, 2, 2): {"parallel": 73, "chair": 27},
    (2, 3, 2): {"chair": 100},
    (2, 4, 2): {"basket": 100},
    (3, 1, 3): {"parallel": 100},
    (3, 2, 3): {"hybrid3": 100},
    (3, 3, 3): {"hybrid1": 50, "hybrid3": 27, "chair": 19, "basket2": 4},
    (3, 4, 3): {"basket": 75, "basket2": 25},
    (4, 1, 4): {"basket": 14},
    (4, 4, 4): {"basket": 43, "chair": 14},
}


def _table3_lookup(a, b, c):
    """Return the {class: percentage} dict from TABLE3, or None."""
    if a != c:
        return None
    key = (a, b, c)
    if key in TABLE3:
        return TABLE3[key]
    key_any = (a, "any", c)
    if key_any in TABLE3:
        return TABLE3[key_any]
    return None


def viability(l1, l2, l3, topology, n_tetrads):
    """Score whether a query fold is viable per gkag435 topological rules.

    Implements hard geometric constraints (loop type vs length, stem height),
    the two-1-nt-loops rule, Table 3 lookup (L1=L3 cases), two-tetrad paper
    rules as fallback, tetrad-count adjustments, and a general propeller
    fallback. Returns one of ``viable``/``marginal``/``not_viable``/
    ``unknown``.
    """
    if n_tetrads < 2:
        return "n/a"

    a, b, c = _clamp(l1), _clamp(l2), _clamp(l3)
    types = loop_type_pattern(topology)
    classes = classify_topology(topology)

    # --- Hard constraints (geometric impossibilities) ---
    # 0-1 nt loops must be propeller
    for length, loop_type in zip((l1, l2, l3), types):
        if length <= 1 and loop_type != "p":
            return "not_viable"
    # Diagonal loops require >= 4 nts
    for length, loop_type in zip((l1, l2, l3), types):
        if loop_type == "d" and length < 4:
            return "not_viable"
    # 4-tetrad stem (~18 Å) too tall for 1-nt propeller loop (~12.5 Å)
    if n_tetrads >= 4 and any(length <= 1 for length in (l1, l2, l3)):
        return "not_viable"

    # --- Two-1-nt-loops rule: topology must be parallel ---
    short_count = sum(1 for length in (l1, l2, l3) if length <= 1)
    if short_count >= 2:
        if "parallel" in classes:
            return "viable"
        return "not_viable"

    # --- d+pd / d-pd explicit rule (basket with two diagonal flanking loops) ---
    if topology in ("d+pd", "d-pd"):
        return "viable" if a == 4 and c == 4 else "not_viable"

    # --- Table 3 lookup (L1 == L3) ---
    entry = _table3_lookup(a, b, c)
    if entry is not None:
        matched = [pct for cls, pct in entry.items() if cls in classes]
        if matched:
            if max(matched) >= 50:
                label = "viable"
            else:
                label = "marginal"
            if label == "viable":
                label = _tetrad_adjustment(classes, n_tetrads, label)
            return label
        return "not_viable"

    # --- Two-tetrad paper rules (fallback for 2-tetrad, not covered by Table 3) ---
    if n_tetrads == 2:
        result = _two_tetrad_rules(a, b, c, topology)
        if result is not None:
            return result

    # --- General propeller fallback ---
    if all(t == "p" for t in types):
        if all(length >= 2 for length in (l1, l2, l3)):
            label = "viable"
            label = _tetrad_adjustment(classes, n_tetrads, label)
            return label
        return "not_viable"

    # --- Topology-specific fallback for non-propeller with loops >= 2 ---
    if all(length >= 2 for length in (l1, l2, l3)):
        if classes and "unknown" not in classes:
            label = "marginal"
            label = _tetrad_adjustment(classes, n_tetrads, label)
            return label

    return "unknown"


def _tetrad_adjustment(classes, n_tetrads, label):
    """Downgrade viable to marginal for energetically unfavorable tetrad counts."""
    if label != "viable":
        return label
    if n_tetrads == 3 and classes & {"chair", "basket", "basket2"}:
        return "marginal"
    if n_tetrads == 2 and classes & {"parallel", "hybrid1", "hybrid2", "hybrid3"}:
        return "marginal"
    return label


def _two_tetrad_rules(a, b, c, topology):
    """Two-tetrad-specific rules from the companion paper. Returns label or None.

    Only handles topologies with explicit rules not well-covered by Table 3:
    -ld+l, +ld-l, -l-l-l. The d+pd case is handled explicitly before Table 3,
    and +l+l+l (chair) is left to Table 3.
    """
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

    return None


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
