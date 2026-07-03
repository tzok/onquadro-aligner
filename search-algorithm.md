# The matching search algorithm

This document explains how `02-find.py` matches a query sequence against a
template quadruplex pattern, why the naive approach is combinatorially
explosive, and how the default search reformulates the problem as a dynamic
program that returns **exactly the same result**, but in polynomial instead
of combinatorial time.

## 1. The problem

Each entry in `onquadro-aligner.json` (a "template") describes one
deduplicated G-quadruplex pattern observed in the PDB, in the form of:

- `qrs`: a list of `4 * n_tetrads` integers, one per G-tetrad column, telling
  which of the `n_tetrads` tetrads each of the template's guanines belongs to
  (only used to reconstruct display/viability metadata — see below).
- `description`: a tuple of `4 * n_tetrads + 1` strings — the nucleotide
  stretches that sit **before**, **between**, and **after** the template's
  own tetrad guanines. These are the "tracts"/"loops" that a matching query
  sequence is expected to reproduce.

Given a query `sequence` and a desired tetrad count `n`, the tool must choose
`k = 4 * n` of the query's G positions (in increasing order) as the
"tetrad-forming" guanines, then compare the `k + 1` segments of the query
that fall *between* consecutive chosen Gs (and before the first / after the
last) against the template's `k + 1` `description` segments. The chosen
combination that minimizes the mismatch is the best match for that template.

The mismatch is `(score_tract, linker_distance)`, compared **lexicographically**
(`score_tract` first, `linker_distance` as a tie-breaker) — see
`process_single_item`/`match_quadruplexes` in `02-find.py`.

## 2. The brute-force approach (`--brute-force`)

The original implementation (still available via `--brute-force`, in
`describe()` + `score_description()`) does exactly what the problem
statement says, literally:

```python
for combination in itertools.combinations(g_indices, k):
    description = describe(sequence, obj["qrs"], combination)     # split query at `combination`
    score_tract, linker_distance, _ = score_description(obj["description"], description)
    ...
```

For a query with `m` G's and `k = 4 * n_tetrads` required tetrad positions,
this enumerates `C(m, k)` combinations — **for every template with that
tetrad count**. For the example that motivated this work (35-nt query, 18
Gs, 3 tetrads → `k=12`), `C(18, 12) = 18,564` combinations, repeated across
**201** three-tetrad templates: ~3.7 million combinations, each running up
to 4 Biopython pairwise alignments. Measured: **~5.15 seconds per template**,
**~110 seconds total** even parallelized across 12 CPU cores. The cost grows
combinatorially with the number of Gs in the query — a longer query, or one
with more spurious (non-quadruplex) Gs, quickly becomes intractable.

## 3. The key insight: the score decomposes over segments

Look closely at `describe()`: it marks the chosen combination indices in the
query, then splits the sequence into segments **wherever a chosen index
is**, regardless of *which* tetrad code was assigned to it:

```python
qrs = [0] * len(sequence)
for code, i in zip(quadruplex_qrs, combination):
    qrs[i] = code          # code is only used to mark "cut here"; its value is discarded below
...
for i, _ in enumerate(qrs):
    if qrs[i] == 0:
        current.append(sequence[i])
    else:
        description.append(current.copy())   # segment boundary — code value irrelevant
        current = []
```

In other words: **the split of the query into segments depends only on
*which* G positions were chosen as cut points, not on the per-tetrad
labelling**. And `score_description()` sums independent, per-segment
contributions:

- `score_tract`: for each *interior* segment (not the first or last), if the
  **template's** segment at that position is empty, add the length of the
  **query's** segment there (a "tract" that should have been a direct
  G-stack got extra nucleotides inserted).
- `linker_distance`: for each segment (including first/last), if both sides
  are non-empty, run a global pairwise alignment and add its (negated)
  score; if only one side is empty, add the other side's length; if both
  are empty, add 0.

Crucially, **the cost of segment `i` depends only on the substring of the
query between cut point `i-1` and cut point `i`, and on the template's fixed
segment `i`** — it has no interaction with any other segment's cost. This
is exactly the structure of an **optimal segmentation / shortest-path
problem**, and it is what makes dynamic programming applicable: we don't
need to enumerate whole combinations, just find, incrementally, the best way
to reach each possible "last cut point" after placing `j` cuts.

## 4. The dynamic program (`find_best_combination_dp`)

### State

Let `positions = [-1, g_0, g_1, ..., g_{m-1}]` be the query's G indices, with
a prepended virtual position `-1` representing "before the sequence starts"
(`positions_with_virtual_start`). Define:

```
dp[j][idx] = the lexicographically-smallest (score_tract, linker_distance)
             achievable by placing the first j cut points, with the j-th
             (last) cut point at positions[idx]
```

`idx` ranges over `0..m` (`0` is the virtual start, only valid for `j=0`).

### Transition

To place the `j`-th cut point at `positions[idx_b]`, we must have come from
some earlier cut point `positions[idx_a]` with `idx_a < idx_b` (or the
virtual start, if `j=1`). The segment consumed by this transition is
`sequence[positions[idx_a]+1 : positions[idx_b]]` — the query nucleotides
strictly between the two chosen Gs (the Gs themselves are excluded, matching
`describe()`). We compare that segment against `template[j-1]` (the
`j`-th template segment, 0-indexed) using `segment_cost()`:

```
dp[j][idx_b] = min over valid idx_a of:
    dp[j-1][idx_a] ⊕ segment_cost(sequence[pos_a+1:pos_b], template[j-1], interior)
```

where `⊕` is component-wise addition of the `(tract, linker)` pair, `min` is
the lexicographic minimum, and `interior = (1 <= j-1 <= k-1)` matches the
`score_tract` rule (first/last segments never contribute to tract distance).

### Base case and final answer

`dp[0][0] = (0, 0)` (no cuts placed yet, zero cost); every other `dp[0][idx]`
is `(∞, ∞)` (you cannot have already placed a cut anywhere without spending a
transition). After filling `dp[1..k]`, the **last** segment (from the k-th
cut point to the end of the sequence) is not yet accounted for — it is added
in one final step outside the main loop, comparing `sequence[pos_b+1:]`
against `template[k]` (the last template segment):

```
answer = min over idx_b of dp[k][idx_b] ⊕ segment_cost(sequence[pos_b+1:], template[k], interior=False)
```

### Reconstructing the combination

Because the search needs to return an actual combination of query indices
(for QRS rendering, viability assessment, and g4composer export — not just a
score), each `dp[j][idx_b]` transition records its best predecessor in a
`back[j][idx_b]` table. After finding the best final `idx_b`, we walk `back`
from `j=k` down to `j=1` to recover the `k` chosen positions, then map them
back through `positions` to get the actual query indices — this **is** the
`combination` tuple that the rest of the pipeline (QRS lines, viability,
g4composer) already expects, byte-for-byte identical in shape to what
`itertools.combinations` would have produced.

### Complexity

- **States**: `O(k * m)` (`k+1` cut-count values times `m+1` positions).
- **Transitions per state**: `O(m)` (any earlier position is a candidate
  predecessor).
- **Total**: `O(k * m²)` segment-cost evaluations, each doing at most one
  Biopython pairwise alignment (itself cached by `(template_segment,
  candidate_segment)` via `functools.cache` on `segment_cost`/`align`, so
  repeated substrings across different `j` are only aligned once).

Compare to brute force's `O(C(m, k))`, which is combinatorial in `m`. For
`m=18, k=12`: DP does at most `13 * 18² ≈ 4,200` segment evaluations versus
`18,564` full combinations (each needing up to 4 alignments — ~74,000
potential alignments). The real payoff shows up for larger `m` (queries with
many Gs): e.g. `m=30, k=16` gives `C(30,16) ≈ 145,000,000` brute-force
combinations, versus DP's `17 * 30² ≈ 15,300` — the DP's advantage grows
without bound as the query gets more Gs, while brute force gets
combinatorially worse.

### Correctness: why lexicographic-min DP is valid here

Dynamic programming with a lexicographically-ordered, vector-valued cost is
only correct if the ordering is compatible with the additive transitions —
i.e. Bellman's principle of optimality must hold for the *lexicographic*
comparison, not just a scalar sum. This holds here because `score_tract` and
`linker_distance` are **summed independently** (never traded off against
each other within a single segment's cost) and are both non-negative:

```
if (a1, b1) < (a2, b2) lexicographically, then for any (c, d):
    (a1 + c, b1 + d) < (a2 + c, b2 + d) lexicographically
```

This is straightforward: if `a1 < a2`, then `a1 + c < a2 + c` regardless of
`b1, b2, d`. If `a1 == a2` and `b1 <= b2`, then `a1 + c == a2 + c` and
`b1 + d <= b2 + d`. So taking the lexicographically-smallest partial cost at
every state and extending it is guaranteed to produce the
lexicographically-smallest total cost — exactly what
`match_quadruplexes`'s final `(score_tract, linker_distance) < (...)`
comparison expects.

This was also verified empirically: the DP path was cross-checked against
brute force on multiple templates (varying tetrad counts and score
profiles), producing byte-identical `(score_tract, linker_distance,
combination)` results in every case, and a full regression across all 201
three-tetrad templates in the database showed zero differences.

## 5. Searching every tetrad count by default

Before the DP search existed, `main()` defaulted to a single hypothesis:
`tetrad_count = [max_tetrad_count]`, where `max_tetrad_count = len(g_indices)
// 4` — the largest number of tetrads the query's *total* G count could
possibly form. This was a deliberate trade-off to keep the brute-force
search's combinatorial cost down (fewer, larger tetrad counts mean *fewer*
Gs are left unused as cut points, which shrinks `C(m, k)`) — but it meant
that if the correct template actually had **fewer** tetrads than
`max_tetrad_count` (i.e. the query contains extra, non-quadruplex G's), it
was never even attempted. This was exactly the original bug report: a
35-nt, 18-G query (`max_tetrad_count = 4`) whose true match was a 3-tetrad
template, which the default search silently skipped.

Now that the DP search makes each tetrad-count hypothesis cheap regardless
of how many spurious Gs the query has, `main()` instead defaults to
searching **every tetrad count the database contains**, from the smallest
up to `max_tetrad_count`:

```python
min_tetrad_count = min((len(d["qrs"]) // 4 for d in data), default=max_tetrad_count)
if max_tetrad_count >= min_tetrad_count:
    tetrad_count = list(range(min_tetrad_count, max_tetrad_count + 1))
else:
    tetrad_count = [max_tetrad_count]
```

`process_single_item`/`match_quadruplexes` already accepted `tetrad_count`
as a list and iterated over it (each database entry only ever matches its
own actual tetrad count — see the `if len(obj["qrs"]) // 4 != n: continue`
guard), so no change was needed there: sweeping multiple tetrad counts
simply means more entries produce a candidate row instead of being skipped.
The final ranking (perfect match first, then viability, then tract
distance, then linker distance — in `main()`'s `result.sort(...)` call)
already compares candidates from different tetrad counts on equal footing,
so a 3-tetrad perfect match correctly outranks a 4-tetrad imperfect one
without any extra logic.

`--tetrad-count N` still restricts the search to exactly one hypothesis,
useful when you already know the expected fold and want to skip the (still
cheap, but not free) work of trying the others.

## 6. Edge cases

- **Fewer available Gs than required** (`k > m`): `find_best_combination_dp`
  returns `None` immediately, matching brute force's behavior of producing
  no combinations at all (`itertools.combinations` of a sample smaller than
  the requested size yields nothing).
- **Unreachable states**: `dp` cells that are still `(∞, ∞)` (no valid
  predecessor existed, e.g. because there weren't enough earlier Gs to place
  `j-1` prior cuts) are skipped when scanning for predecessors/final answers,
  so `∞` never leaks into a reported score.
- **Query too short for any database entry** (`max_tetrad_count < min_tetrad_count`):
  falls back to `tetrad_count = [max_tetrad_count]` (an empty/degenerate
  search that correctly reports no matches) rather than an empty `range()`.

## 7. Where to look in the code

| Concept                                    | Function                            |
| ------------------------------------------- | ------------------------------------ |
| Per-segment cost (tract + linker)           | `segment_cost()`                    |
| The DP search itself                        | `find_best_combination_dp()`        |
| Virtual "before sequence start" position    | `positions_with_virtual_start()`    |
| Brute-force equivalent (kept for debugging) | `describe()` + `score_description()`|
| Dispatch between DP and brute force         | `process_single_item()`             |
| Default tetrad-count sweep                  | `main()`                            |


## 7. When to still use `--brute-force`

The DP and brute-force paths are verified to produce identical results, so
there is normally no reason to prefer brute force. It is kept as a CLI flag
(`--brute-force`) purely as a slower, more "obviously correct by inspection"
reference implementation — useful when validating future changes to the
scoring rules, or when debugging a suspected discrepancy between the two
search strategies.
