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

**QRS notation**: Letters encode G-tetrad stacking. The same letter marks all four Gs of one tetrad. Dots (`.`) are non-tetrad nucleotides.

### Options

| Flag                         | Description                                              |
| ---------------------------- | -------------------------------------------------------- |
| `--molecule DNA\|RNA\|Other` | Filter by molecule type                                  |
| `--tetrad-count N`           | Look for N tetrads instead of the maximum possible       |
| `--list-all-quadruplex`      | Show all matching file names instead of only one per hit |
| `--output-csv FILE`          | Write detailed results (incl. viability) to a CSV file  |
| `--g4composer-output-dir DIR`| Write g4composer `.inp` files (one per matched template) |
| `--jobs N`                   | Number of parallel workers (default: all CPUs)           |
| `--input FILE`               | Use a different preprocessed JSON database               |

### g4composer export

```bash
uv run ./02-find.py GGTTGGCGCGAAGCATTCGCGGGTTGG --tetrad-count 2 --g4composer-output-dir ./g4c/
```

Writes one `.inp` file per matched template, named `{rank}-{template}.inp` where `{rank}` is the zero-padded position of the match in the sorted result list (best match = `0`), so files sort by overall match quality. The sort key is `(perfect, viability, tract distance, linker distance)`, so `{rank}` encodes all components in a single number. Each file is a valid g4composer input describing the **query sequence** folded into the matched template's geometry:

``` text
# tract_distance=0 linker_distance=0 viability=unknown loop_lengths=6-1-3 topology=d-p-l template=5zev-assembly1
name        5zev-assembly1
sequence    gggtacccgggtgaggtgcggggt
structure   ^^......^^^.^.^^...^^^^.
chi         ........................
sugar       ........................
orient      A-;B+;C-
rise        -2.9;6.3
twist       -43.9;75.6
path        A1;B1;B3;A3;C3;B2;A2;C2;C1;B4;A4;C4
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
| `n/a`        | Not assessable (no loop info, or <2 tetrads)                            |

**Caveats**: the rules are derived from structural studies of DNA G4s, mostly with thymidine-only loops and no flanking residues. They are applied to query loop **lengths** regardless of composition — non-T loops and flanking sequences can shift conformational preferences (e.g. via base-pair or base-triple formation within loops). Parallel structures resolved in PEG/crowding conditions were excluded from Table 3, so parallel topologies may be underrepresented.

### Rebuild the database

```bash
uv run ./01-preprocess.py /path/to/eltetrado/json/ --output onquadro-aligner.json
```

The input directory should contain ElTetrado JSON output files (one per PDB assembly). The output is a JSON file where each row describes a unique single-chain G-only quadruplex pattern.

## How it works

1. **Preprocessing** (`01-preprocess.py`) — Reads ElTetrado JSON files, extracts single-chain all-G quadruplexes, and deduplicates them into a compact pattern database. Each entry stores the QRS codes, QRS characters, linker descriptions (tract sequences between tetrads), source filenames, a `g4c` field with g4composer geometry (`orient`, `rise`, `twist`, `path`) reconstructed from the ElTetrado JSON topology, and a `loops` field with loop topology and gap indices for viability scoring.

2. **Search** (`02-find.py`) — Takes a nucleotide sequence, identifies G positions, and tries all combinations of Gs to match against database patterns. Each match is scored for topological viability from the query's inferred loop lengths and the template's loop topology. Two distance metrics are used: **tract distance** (extra nucleotides in tracts) and **linker distance** (edit distance between linker sequences, penalizing mismatches, gaps, and missing/extra linkers). Perfect matches (both distances = 0) are ranked first; all others are sorted by viability, then by tract distance, then by linker distance. Searches in parallel across all CPU cores. When `--g4composer-output-dir` is set, writes one `.inp` file per matched template using the query sequence and the template's geometry.

## Database

The pre-built `onquadro-aligner.json` contains **327** unique single-chain G-quadruplex patterns extracted from the ElTetrado corpus. Includes only quadruplexes that:

- Span a single chain (no inter-chain quadruplexes)
- Contain exclusively G-tetrads
