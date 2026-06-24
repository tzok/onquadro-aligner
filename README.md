# onquadro-aligner

Search an RNA or DNA sequence against a database of known single-chain G-quadruplex structures.

The tool matches the input sequence against preprocessed ElTetrado data and outputs **QRS** lines — a dot-bracket-like notation showing the positions and orientations of G-tetrads. It can also export **g4composer** input files for 3D structure modeling.

## Setup

Requires Python 3.14+.

```bash
uv sync
```

Dependencies: `biopython`, `pandas`, `tqdm`.

## Usage

### Search a sequence

```bash
uv run ./02-find.py GGGTACCCGGGTGAGGTGCGGGGT
```

Output:

```
GGGTACCCGGGTGAGGTGCGGGGT
qR......RqS.r.Qs...SrQs.
```

The first line is the input sequence. Each subsequent line is a de-duplicated QRS representation of a matching quadruplex, ordered by best match first. Results are sorted by **tract distance** (lower is better), then by **linker score** (higher is better), then by **viability** (viable first) as a tie-breaker.

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

Writes one `.inp` file per matched template, named `{rank}-{template}.inp` where `{rank}` is the zero-padded position of the match in the sorted result list (best match = `0`), so files sort by overall match quality. The sort tuple is `(tract distance, -linker score, viability)`, so `{rank}` encodes all three in a single number. Each file is a valid g4composer input describing the **query sequence** folded into the matched template's geometry:

``` text
# tract_distance=0 linker_score=38 viability=unknown loop_lengths=2-15-2 topology=+l-l+l template=6fc9-assembly-1
name        6fc9-assembly-1
sequence    ggttggcgcgaagcattcgcgggttgg
structure   ^^..^^...............^^..^^
chi         ...........................
sugar       ...........................
orient      A+;B-
rise        3.4
twist       18.7
path        A1;B1;B2;A2;A3;B3;B4;A4
```

The `# score` line prepends the aligner's tract distance, linker score, and **viability** assessment. `sequence` and `structure` describe the query (DNA lowercase, RNA uppercase; `^` marks the matched tetrad Gs). `orient`, `rise`, `twist`, and `path` come from the reference template and are reconstructed from the ElTetrado JSON during preprocessing — no 3D coordinates required.

### Viability scoring

Each match is assessed for **topological viability** based on the relationship between loop lengths and G4 topology described in a comprehensive study of two-tetrad DNA G-quadruplexes with thymidine-only loops. The four topologies characterized in that work (`+l+l+l`, `−ld+l`, `−l−l−l`, `d+pd`) are applied to any tetrad count. Propeller-containing topologies not among these four are treated as viable for ≥3 tetrads (propeller loops are abundant in three-tetrad G4s) and not viable for two tetrads. Results are sorted by viability first:

| Label        | Meaning                                                                 |
| ------------ | ----------------------------------------------------------------------- |
| `viable`     | Query loop lengths fall in the template topology's preferred region     |
| `marginal`   | Equilibrium/mix region — unimolecular fold possible but not dominant    |
| `unknown`    | Topology not covered by the manuscript's rules                          |
| `not_viable` | Loop lengths incompatible with the template topology (e.g. single-nt)   |
| `n/a`        | Not assessable (no loop info, or <2 tetrads)                            |

**Caveat**: the rules were derived for two-tetrad DNA G4s with thymidine-only loops and no flanking residues. They are applied to query loop **lengths** regardless of composition — non-T loops and flanking sequences can shift conformational preferences (e.g. via base-pair or base-triple formation within loops).

### Rebuild the database

```bash
uv run ./01-preprocess.py /path/to/eltetrado/json/ --output onquadro-aligner.json
```

The input directory should contain ElTetrado JSON output files (one per PDB assembly). The output is a JSON file where each row describes a unique single-chain G-only quadruplex pattern.

## How it works

1. **Preprocessing** (`01-preprocess.py`) — Reads ElTetrado JSON files, extracts single-chain all-G quadruplexes, and deduplicates them into a compact pattern database. Each entry stores the QRS codes, QRS characters, linker descriptions (tract sequences between tetrads), source filenames, a `g4c` field with g4composer geometry (`orient`, `rise`, `twist`, `path`) reconstructed from the ElTetrado JSON topology, and a `loops` field with loop topology and gap indices for viability scoring.

2. **Search** (`02-find.py`) — Takes a nucleotide sequence, identifies G positions, and tries all combinations of Gs to match against database patterns. Each match is scored for topological viability from the query's inferred loop lengths and the template's loop topology. Ranking is by viability, then by tract distance (lower is better), then by linker alignment score (higher is better). Searches in parallel across all CPU cores. When `--g4composer-output-dir` is set, writes one `.inp` file per matched template using the query sequence and the template's geometry.

## Database

The pre-built `onquadro-aligner.json` contains **327** unique single-chain G-quadruplex patterns extracted from the ElTetrado corpus. Includes only quadruplexes that:

- Span a single chain (no inter-chain quadruplexes)
- Contain exclusively G-tetrads
