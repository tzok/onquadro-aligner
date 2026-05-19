# onquadro-aligner

Search an RNA or DNA sequence against a database of known single-chain G-quadruplex structures.

The tool matches the input sequence against preprocessed ElTetrado data and outputs **QRS** lines — a dot-bracket-like notation showing the positions and orientations of G-tetrads.

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

The first line is the input sequence. Each subsequent line is a de-duplicated QRS representation of a matching quadruplex, ordered by best match first.

**QRS notation**: Letters encode G-tetrad stacking. The same letter marks all four Gs of one tetrad. Dots (`.`) are non-tetrad nucleotides.

### Options

| Flag                         | Description                                              |
| ---------------------------- | -------------------------------------------------------- |
| `--molecule DNA\|RNA\|Other` | Filter by molecule type                                  |
| `--tetrad-count N`           | Look for N tetrads instead of the maximum possible       |
| `--list-all-quadruplex`      | Show all matching file names instead of only one per hit |
| `--output-csv FILE`          | Write detailed results to a CSV file                     |
| `--jobs N`                   | Number of parallel workers (default: all CPUs)           |
| `--input FILE`               | Use a different preprocessed JSON database               |

### Rebuild the database

```bash
uv run ./01-preprocess.py /path/to/eltetrado/json/ --output onquadro-aligner.json
```

The input directory should contain ElTetrado JSON output files (one per PDB assembly). The output is a JSON file where each row describes a unique single-chain G-only quadruplex pattern.

## How it works

1. **Preprocessing** (`01-preprocess.py`) — Reads ElTetrado JSON files, extracts single-chain all-G quadruplexes, and deduplicates them into a compact pattern database. Each entry stores the QRS codes, QRS characters, linker descriptions (tract sequences between tetrads), and source filenames.

2. **Search** (`02-find.py`) — Takes a nucleotide sequence, identifies G positions, and tries all combinations of Gs to match against database patterns. Ranking is by tract distance (lower is better), then by linker alignment score (higher is better). Searches in parallel across all CPU cores.

## Database

The pre-built `onquadro-aligner.json` contains **327** unique single-chain G-quadruplex patterns extracted from the ElTetrado corpus. Includes only quadruplexes that:

- Span a single chain (no inter-chain quadruplexes)
- Contain exclusively G-tetrads
