# The Secret Trick: Why △ is Deterministic

The entire computational engine of `lr_coefficients.py` rests on one
structural observation about KTW puzzle tiles that isn't obvious from
the paper or from the brute-force solver.

## The observation

There are two kinds of unit triangles in a puzzle: upward (△) and downward (▽).

**△ is deterministic.**  Given the bottom edge label and the left╱ edge label,
the right╲ edge label is *uniquely forced* by the Delta triple constraints.

**▽ is where all the branching happens.**  Given the left╲ edge label (which
came from the △ to its left), there are *at most 2* valid choices of
(top edge, right╱ edge).

This means:

| left╲ | choices for (top, right╱) |
|-------|---------------------------|
| 0     | (0, 0) or (10, 1)         |
| 1     | (0, 10) or (1, 1)         |
| 10    | (1, 0) — only one!        |

## Why this matters

The brute-force solver (`solve_puzzle.py`) treats all internal edges as
unknowns and tries all 3^n combinations.  For a side-N triangle, that's
exponential in N².

The transfer matrix approach exploits the determinism of △:

1. Process each row left-to-right
2. At each △: the right╲ diagonal is forced (one lookup, no branching)
3. At each ▽: branch into ≤2 possibilities
4. The *effective* branching factor per row is at most 2^(w-1), where w
   is the row width — and in practice much less, because many branches
   die at the next △ (invalid Delta pair)

The row transfer matrix is thus a *linear-time scan* per row with
bounded branching, rather than an exponential enumeration.

## The connection to commutativity

The row transfer decomposes into two commuting pieces:
- **T̃_+**: processes the right boundary (red spectral line)
- **T̃_-**: processes the left boundary (green spectral line)

These commute because of the Yang-Baxter equation for the free-fermion
vertex model (Zinn-Justin 2008, Theorem 4).  In practice, this means:

- You can process left boundary edges and right boundary edges in *any order*
- The puzzle count (= LR coefficient) is independent of the processing order
- Each processing order gives a different *natural basis* for the puzzle space,
  but they're all in bijection

The row-by-row approach interleaves one T̃_+ and one T̃_- per row.
But you could equally apply all T̃_+ first, then all T̃_- — and get the
same count.

## The △ determinism table

```
Delta triples: (bottom, left╱) → right╲
  (0, 0)  →  0
  (1, 0)  →  10
  (10, 1) →  0
  (1, 1)  →  1
  (0, 10) →  1

Invalid pairs (no valid right╲):
  (0, 1), (1, 10), (10, 0), (10, 10)
```

An invalid pair at △ kills the entire branch — this is how the transfer
matrix prunes the search space aggressively.

## The ▽ branching table

```
Nabla: given left╲, valid (top, right╱) pairs:
  left╲ = 0  →  (0, 0)  or  (10, 1)      [2 choices]
  left╲ = 1  →  (0, 10) or  (1, 1)       [2 choices]
  left╲ = 10 →  (1, 0)                    [1 choice — forced!]
```

Note: the rhombus label (10) entering a ▽ forces a unique outcome,
just like △.  So branches that produce many rhombi converge quickly.

## Reconstruction

Because △ is deterministic and ▽ is determined given its top edge label
(from the next row's horizontal edges), the *horizontal edge labels at
each level* are a complete description of the puzzle.  All diagonal labels
can be reconstructed in O(N²) time.  This is what `Puzzle.edge_labels()` does.
