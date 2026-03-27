# KTW Puzzles & Littlewood-Richardson Coefficients

Three approaches to computing Littlewood-Richardson coefficients c^ν_{λ,μ} via Knutson-Tao-Woodward puzzle triangles, based on Zinn-Justin (2008) ["Littlewood-Richardson coefficients and integrable tilings"](https://arxiv.org/abs/0809.2392).

## Quick start

```python
from lr_coefficients import lr_count, lr_puzzles

# The famous c = 2 case
lr_count((2,1), (2,1), (3,2,1))   # → 2
lr_puzzles((2,1), (2,1), (3,2,1)) # → [Puzzle, Puzzle]

# Hopf algebra coproduct
from transfer_operators import lr_via_transfer

coprod = lr_via_transfer((3,2,1))  # → {((2,1),(2,1)): 2, ...} (25 terms)
```

Run demos: `python3 lr_coefficients.py` or `python3 transfer_operators.py`

## What's here

### `lr_coefficients.py` — Puzzle enumeration

Row transfer matrix that exploits a key structural fact: upward triangles (△) are **deterministic** — given (bottom, left╱), the right╲ label is forced. All multiplicity comes from downward triangles (▽), which branch at most 2 ways. This turns exponential brute-force into a bounded-branching scan. See [SECRET_TRICK.md](SECRET_TRICK.md).

### `transfer_operators.py` — Fock space operators & coproduct

Fermionic Fock space operators T_free (removes horizontal strips) and T̄_free (removes vertical strips), implementing the free-fermion vertex model from Zinn-Justin.

Includes `lr_via_transfer(lam)` which computes the full **Hopf algebra coproduct**:

```
Δ(s_λ) = Σ_{μ,ν} c^λ_{μ,ν}  s_μ ⊗ s_ν
```

Uses Kostka matrix inversion (unitriangular in dominance order) over T_free chains to extract individual c^λ_{μ,ν}. Cross-validated against the puzzle solver.

### `solve_puzzle.py` — Brute force

Enumerates all internal edge labelings. Exponential but useful for cross-validation.

### LaTeX

- `puzzles.tex` — Edge-labeled puzzle diagrams
- `lr_paths3.tex` — Lattice path visualizations (green/red paths through tile interiors)

## A discovery along the way

The naive approach — evaluate s_{λ/μ}(1,...,1) at multiple n and solve a linear system — **fails** because hook-content polynomials have genuine linear dependencies:

```
dim_n(3,1) − 3·dim_n(2,2) + dim_n(2,1,1) ≡ 0   for all n
```

The principal specialization (all variables equal) loses information about Schur functions. The Kostka matrix, which tracks SSYT by content rather than just total dimension, has no such degeneracy.

## References

- Knutson, Tao, Woodward — "The honeycomb model of GL_n(C) tensor products II" (2004)
- Zinn-Justin — "Littlewood-Richardson coefficients and integrable tilings" ([arXiv:0809.2392](https://arxiv.org/abs/0809.2392))
- Molev — Lecture slides on puzzle combinatorics (included as `molev.pdf`)
