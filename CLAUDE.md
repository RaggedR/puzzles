# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Research on KTW (Knutson-Tao-Woodward) puzzle triangles and Littlewood-Richardson coefficients, based on Zinn-Justin (2008) "Littlewood-Richardson coefficients and integrable tilings" (arXiv:0809.2392).

Reference paper: `~/git/research/data/puzzles/zinn-justin-2008-lr-coefficients-integrable-tilings.pdf`

## Files

### `transfer_operators.py`

Algebraic approach: transfer matrix operators in the fermionic Fock space.

- **Partition ↔ binary string**: `partition_to_binary`, `binary_to_partition` (convention: 1=green/vertical, 0=red/horizontal)
- **T_free**: removes horizontal strips (interlacing condition)
- **T_bar_free**: removes vertical strips (at most 1 box per row)
- **T̃ operators on G_free**: `T_tilde_plus` (T_free on ν), `T_tilde_minus` (T̄_free on μ)
- **Operator algebra**: `apply_op`, `state_single`, `inner_product` for composing operators on formal linear combinations
- **Coproduct**: `lr_via_transfer(lam)` computes the full Hopf algebra coproduct Δ(s_λ) = Σ c^λ_{μ,ν} s_μ ⊗ s_ν using Kostka matrix inversion over T_free chains. Cross-validated against the puzzle solver in `lr_coefficients.py`.
- **Verification**: `schur_via_transfer`, `hook_content_formula` confirm Schur function identities

Run with: `python3 transfer_operators.py` (runs `demo()` with all verifications).

### `lr_coefficients.py`

Commuting transfer operator approach: computes c^ν_{λ,μ} and enumerates KTW puzzle witnesses via the row transfer matrix on trinary-labeled edges.

- **Row transfer matrix**: processes △ (deterministic) and ▽ (branching) triangles left-to-right
- **`lr_puzzles(lam, mu, nu)`**: returns list of `Puzzle` objects — one per puzzle witness
- **`lr_count(lam, mu, nu)`**: just the coefficient (faster, uses `Counter`)
- **`Puzzle.edge_labels()`**: reconstructs full edge labeling for visualization
- Key structural fact: △ is deterministic given (bottom, left╱); all multiplicity comes from ▽ branching

Run with: `python3 lr_coefficients.py` (runs `demo()` with verification against known LR values).

### `solve_puzzle.py`

Brute-force approach: assembles KTW puzzle triangles by enumerating all internal edge labelings.

### `puzzles.tex` / `puzzles.pdf`

LaTeX writeup with puzzle triangle diagrams.

### Working Notes

- `encoding-notes.md` — Binary encoding of partitions (formula: position of i-th 1 = m + i - λᵢ)
- `boundary-computation.md` — Tracing the geometry of example ZJ puzzles
- `LR-G.pdf` — Reference diagram

## Development Direction

- **Coproduct ↔ puzzles bridge is closed.** `lr_via_transfer` computes c^λ_{μ,ν} purely from the Fock space operators, `lr_count`/`lr_puzzles` from puzzle combinatorics — both give the same answer. The coproduct decomposes s_{(3,2,1)} into 25 tensor terms with the famous c=2 at (μ,ν)=((2,1),(2,1)).
- **Key discovery**: the "evaluate at (1,...,1)" approach fails because hook-content polynomials have linear dependencies — e.g., dim_n(3,1) − 3·dim_n(2,2) + dim_n(2,1,1) ≡ 0. The fix uses Kostka matrix inversion (unitriangular in dominance order).
- Visualization: draw puzzle objects with colored lattice paths (partially done in `lr_paths3.tex`).
- Yang-Baxter commutativity: implement T̃_± as separate operators and verify T̃_+ T̃_- = T̃_- T̃_+ explicitly.
- **Open question**: is there a Hopf algebra of dissections where the stomachion's Sol = ∐ ∏ Sol(Rⱼ) is literally the coproduct formula? See `~/git/wooden/speculative/` for the stomachion project.
