# onquadro-aligner

Search an RNA or DNA sequence against a database of known single-chain G-quadruplex structures.

The tool matches the input sequence against preprocessed ElTetrado data and outputs **QRS** lines — a dot-bracket-like notation showing the positions and orientations of G-tetrads. It can also export **g4composer** input files for 3D structure modeling.

## Setup

Requires Python 3.14+.

```bash
uv sync
```

Dependencies: `biopython`, `pandas`, `streamlit`, `tqdm`.

## Usage

### Web app

```bash
uv run streamlit run app.py
```

Enter a sequence in the input field (uppercase = RNA, lowercase = DNA, mixed = search both). Results appear in a sortable table with viability indicators. Click a row to view its g4composer output and download the `.inp` file.

### Search a sequence (CLI)

```bash
uv run ./02-find.py GGGTACCCGGGTGAGGTGCGGGGT
```

Output:

```
GGGTACCCGGGTGAGGTGCGGGGT
qR......RqS.r.Qs...SrQs.
```

The first line is the input sequence. Each subsequent line is a de-duplicated QRS representation of a matching quadruplex, ordered by best match first. **Perfect matches** (tract distance = 0 and linker distance = 0) are always ranked first. All other matches are sorted by **viability** (viable first), then by **tract distance** (lower is better), then by **linker distance** (lower is better).

By default, the search tries **every tetrad count present in the database** (not just the maximum the query's G count could support), so a template with fewer tetrads than the query has G's for can still win — e.g. a 35-nt query with 18 G's (which could in principle form 4 tetrads) is also matched against 2- and 3-tetrad templates, and the 3-tetrad template turns out to be the correct, perfect match. Use `--tetrad-count` to restrict the search to one specific count.

**QRS notation**: Letters encode G-tetrad stacking. The same letter marks all four Gs of one tetrad. Dots (`.`) are non-tetrad nucleotides.

### Options

| Flag                         | Description                                              |
| ---------------------------- | -------------------------------------------------------- |
| `--molecule DNA\|RNA\|Other` | Filter by molecule type                                  |
| `--tetrad-count N`           | Look for exactly N tetrads instead of searching every tetrad count the database supports |
| `--list-all-quadruplex`      | Show all matching file names instead of only one per hit |
| `--output-csv FILE`          | Write detailed results (incl. viability) to a CSV file  |
| `--g4composer-output-dir DIR`| Write g4composer `.inp` files (one per matched template) |
| `--jobs N`                   | Number of parallel workers (default: all CPUs)           |
| `--input FILE`               | Use a different preprocessed JSON database               |
| `--brute-force`              | Use the exhaustive combinatorial search instead of the default dynamic-programming search (see [search-algorithm.md](search-algorithm.md)); much slower, kept for debugging/verification |

### g4composer export

```bash
uv run ./02-find.py GGTTGGCGCGAAGCATTCGCGGGTTGG --tetrad-count 2 --g4composer-output-dir ./g4c/
```

Writes one `.inp` file per matched template, named `{rank}-{template}.inp` where `{rank}` is the zero-padded position of the match in the sorted result list (best match = `0`), so files sort by overall match quality. The sort key is `(perfect, viability, tract distance, linker distance)`, so `{rank}` encodes all components in a single number. Each file is a valid g4composer input describing the **query sequence** folded into the matched template's geometry:

``` text
# tract_distance=0 linker_distance=0 viability=not_viable loop_lengths=2-15-2 topology=+l+l+l template=2m8z-assembly1
name        2m8z-assembly1
sequence    ggttggcgcgaagcattcgcgggttgg
structure   ^^..^^...............^^..^^
chi         ...........................
sugar       ...........................
orient      A+;B-
rise        3.9
twist       19.2
path        A1;B1;B2;A2;A3;B3;B4;A4
```

The `# score` line prepends the aligner's tract distance, linker score, and **viability** assessment. `sequence` and `structure` describe the query (DNA lowercase, RNA uppercase; `^` marks the matched tetrad Gs). `orient`, `rise`, `twist`, and `path` come from the reference template and are reconstructed from the ElTetrado JSON during preprocessing — no 3D coordinates required.

### Viability scoring

Each match is assessed for **topological viability** based on rules from two studies of intramolecular G-quadruplex folding:

1. **Hard geometric constraints** (gkag435): 0–1 nt loops must be propeller; diagonal loops require ≥4 nts; 4-tetrad stems need ≥2-nt propeller loops.
2. **Two-1-nt-loops rule** (gkag435): when two loops are 1-nt long, the topology must be parallel (all propeller).
3. **Table 3 lookup** (gkag435): when L1=L3, observed topology percentages determine viability (≥50% → viable, 4–49% → marginal, 0% → not_viable).
4. **Two-tetrad rules** (companion paper): explicit rules for `d+pd`, `−ld+l`, `+ld−l`, `−l−l−l` topologies as a fallback for 2-tetrad cases not covered by Table 3.
5. **Tetrad-count adjustments**: 3-tetrad antiparallel and 2-tetrad parallel/hybrid are downgraded to marginal (energetically less favorable gc successions or insufficient stacking).

| Label        | Meaning                                                                 |
| ------------ | ----------------------------------------------------------------------- |
| `viable`     | Query loop lengths fall in the template topology's preferred region     |
| `marginal`   | Equilibrium/mix region — unimolecular fold possible but not dominant    |
| `unknown`    | Topology not covered by the rules                                       |
| `not_viable` | Loop lengths incompatible with the template topology (geometric/rules)  |
| `n/a`        | Not assessable (<2 tetrads, or loop gaps could not be determined)       |

**Caveats**: the rules are derived from structural studies of DNA G4s, mostly with thymidine-only loops and no flanking residues. They are applied to query loop **lengths** regardless of composition — non-T loops and flanking sequences can shift conformational preferences (e.g. via base-pair or base-triple formation within loops). Parallel structures resolved in PEG/crowding conditions were excluded from Table 3, so parallel topologies may be underrepresented.

**Loop positions and lengths**: loop positions come from the template's g4composer `path`, which pins every tetrad G to a tract (column) and a tetrad (letter). Two consecutive template Gs are chain-continuous (stacked) when they share the tract and occupy adjacent tetrad layers; every other connector is a chain discontinuity — a loop between tracts, or a re-entry into the same tract (snapback-like). Non-canonical folds (the chain visiting a tract more than once) therefore have **more than three loops**, and `loop_lengths` lists one length per connector, e.g. `2-2-1-0` for a four-loop fold. For templates whose path is unavailable (see below), loop gaps fall back to ElTetrado annotations and finally to a heuristic inference from the QRS pattern.

**Topology derivation**: for templates without ElTetrado loop-type annotations, the topology (propeller/lateral/diagonal per loop) is derived from the path by mirroring ElTetrado's own classification: opposite tracts → diagonal, flanking nucleotides in the same tetrad → lateral, otherwise propeller. The +/− sign is a pure function of the columns connected (ElTetrado anchors column 1 in 3D and derives the rest around the helix axis): `+` when the 3'-flanking tract is the clockwise-next of the 5'-flanking one, `−` when it is the anticlockwise-next. Path-derived signs agree with the ElTetrado annotations on the whole corpus. Sign subtypes never observed in the ASC-G4 scheme (e.g. `9b`, `7b`, `5b`, `4a` — four corpus entries) carry no common topology name.

**0-nt loops**: a 0-nt connector (two tetrad Gs directly adjacent) is structurally exceptional where the template itself requires nucleotides there — the query then failed to supply a loop the fold needs. Such matches are always marked `not_viable` and penalized in the ranking (+1000 per violation on linker distance), even if they would otherwise be a "perfect" match. A template may also contain *native* 0-nt connectors (snapbacks — a chain re-entering its own tract); reproducing those faithfully with 0 nt is allowed — the length still shows up in `loop_lengths` (e.g. `6-1-3-0`), but without a penalty.

### Rebuild the database

```bash
uv run ./01-preprocess.py /path/to/eltetrado/json/ --output onquadro-aligner.json
```

The input directory should contain ElTetrado JSON output files (one per PDB assembly). The output is a JSON file where each row describes a unique single-chain G-only quadruplex pattern.

## How it works

1. **Preprocessing** (`01-preprocess.py`) — Reads ElTetrado JSON files, extracts single-chain all-G quadruplexes, and deduplicates them into a compact pattern database. Each entry stores the QRS codes, QRS characters, linker descriptions (tract sequences between tetrads), source filenames, a `g4c` field with g4composer geometry (`orient`, `rise`, `twist`, `path`) read from the paired `.inp` file (attached only when the `.inp` path is consistent with the entry's own QRS — assemblies containing several single-chain quadruplexes share one `.inp`), and a `loops` field with the ElTetrado loop topology and gap indices (only for canonical 3-loop folds).

2. **Search** (`02-find.py`) — Takes a nucleotide sequence, identifies G positions, and searches for the best-scoring combination of Gs to match against each database pattern, trying every tetrad count the database supports by default (not just the maximum the query's G count could support) so that fewer-tetrad templates aren't missed. By default this uses a dynamic-programming search that is mathematically equivalent to trying every combination but runs in polynomial instead of combinatorial time (see [search-algorithm.md](search-algorithm.md) for a detailed explanation); the exhaustive combinatorial search is still available via `--brute-force` for verification. Each match is scored for topological viability from the query's loop lengths at every chain-discontinuity connector of the template, together with the template's loop topology (annotated by ElTetrado or derived from the path). Two distance metrics are used: **tract distance** (extra nucleotides in tracts) and **linker distance** (edit distance between linker sequences, penalizing mismatches, gaps, and missing/extra linkers). Perfect matches (both distances = 0) are ranked first; all others are sorted by viability, then by tract distance, then by linker distance. Searches in parallel across all CPU cores. When `--g4composer-output-dir` is set, writes one `.inp` file per matched template using the query sequence and the template's geometry.

## Database

The pre-built `onquadro-aligner.json` contains **340** unique single-chain G-quadruplex patterns extracted from the ElTetrado corpus. Includes only quadruplexes that:

- Span a single chain (no inter-chain quadruplexes)
- Contain exclusively G-tetrads
