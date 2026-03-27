# What We Built (2026-03-20)

## LR Coefficients via Commuting Transfer Operators

### The computation: `lr_coefficients.py`

Computes Littlewood-Richardson coefficients c^ν_{λ,μ} and enumerates
individual KTW puzzle witnesses using the row transfer matrix.

**API:**
```python
from lr_coefficients import lr_puzzles, lr_count

lr_count((2,1), (2,1), (3,2,1))        # → 2
puzzles = lr_puzzles((2,1), (2,1), (3,2,1))  # → [Puzzle, Puzzle]
puzzles[0].levels      # trinary string at each horizontal level
puzzles[0].edge_labels()  # full edge labeling for visualization
```

**The secret trick:** △ (upward) tiles are deterministic — given (bottom, left╱),
the right╲ label is uniquely forced by the Delta triple constraints. All puzzle
multiplicity comes from ▽ (downward) tiles, which branch at most 2 ways.
This turns exponential brute-force search into a row transfer matrix with
bounded branching, heavily pruned by deterministic △ kills.

See `SECRET_TRICK.md` for the full writeup.

**Verified against:**
- 16 hand-computed LR coefficients (Pieri, column×column, row×column, c=2)
- Symmetry c^ν_{λ,μ} = c^ν_{μ,λ}
- Full product s_{(2,1)}² = s_{(4,2)} + s_{(4,1,1)} + s_{(3,3)} + 2s_{(3,2,1)} + s_{(3,1,1,1)} + s_{(2,2,2)} + s_{(2,2,1,1)}
- Cross-validation against brute-force `solve_puzzle.py`

### The visualization: lattice path model

Draws puzzles in the Zinn-Justin/Molev lattice path style (see Molev slides pp 51–53):
- Green and red paths go through triangle **interiors**, connecting edge midpoints
- Each of the 10 tile types (5 △ + 5 ▽) has a specific routing rule
- Paths are traced as continuous polylines from boundary to boundary

**Routing rules per tile:**

△ (bottom─, left╱, right╲):
| Triple     | Green              | Red                |
|------------|--------------------|--------------------|
| (0,0,0)    | —                  | left╱ ↔ right╲     |
| (1,0,10)   | bottom─ ↔ right╲   | left╱ ↔ right╲     |
| (10,1,0)   | bottom─ ↔ left╱    | bottom─ ↔ right╲   |
| (1,1,1)    | left╱ ↔ right╲     | —                  |
| (0,10,1)   | left╱ ↔ right╲     | bottom─ ↔ left╱    |

▽ (top─, right╱, left╲):
| Triple     | Green              | Red                |
|------------|--------------------|--------------------|
| (0,0,0)    | —                  | left╲ ↔ right╱     |
| (0,10,1)   | left╲ ↔ right╱     | top─ ↔ right╱      |
| (10,1,0)   | top─ ↔ right╱      | top─ ↔ left╲       |
| (1,0,10)   | top─ ↔ left╲       | right╱ ↔ left╲     |
| (1,1,1)    | left╲ ↔ right╱     | —                  |

The rule: green connects the 2 edges carrying green (label 1 or 10),
red connects the 2 edges carrying red (label 0 or 10). For all-same
tiles (0,0,0) and (1,1,1), only the diagonal pair is connected; the
horizontal edge's strand comes from the adjacent triangle.

**Generated files:**
- `lr_paths3.tex` / `lr_paths3.pdf` — current best rendering
- `lr_paths2.tex` — per-triangle segments (same visual, different code)
- `lr_paths.tex` — earlier attempt with coloured edges (wrong model)
- `lr_puzzles.tex` — edge-label + old path model views

### Still to do

- Smooth the paths (Bezier curves instead of straight midpoint-to-midpoint segments)
- Offset green/red at label-10 crossings so both are visible
- Match the full Molev page 53 aesthetic (paths extending to boundary, clean line weights)
- Add the drawing code as a proper function in `lr_coefficients.py` or a separate module

### File inventory

| File | Purpose |
|------|---------|
| `lr_coefficients.py` | LR coefficients + puzzle enumeration via transfer matrix |
| `transfer_operators.py` | F-space Schur function operators (T_free, T̄_free) |
| `solve_puzzle.py` | Brute-force puzzle solver + TikZ generation |
| `SECRET_TRICK.md` | Why △ is deterministic |
| `CLAUDE.md` | Project documentation for Claude Code |
| `encoding-notes.md` | Binary encoding derivation |
| `boundary-computation.md` | Hand-traced ZJ example |
| `puzzles.tex` | Original puzzle diagrams |
| `molev.pdf` | Reference: Molev lecture slides |
