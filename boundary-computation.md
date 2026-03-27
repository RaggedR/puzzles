# Boundary Computation for ZJ Example

## Setup
λ = (2,1,1), μ = ν = (1,1), n = 3, m = 2, side length N = 5
Formula: position of i-th 1 = m + i - λ_i (1-indexed)

## λ = (2,1,1) — bottom edge (read left to right)
- i=1: 2+1-2 = 1
- i=2: 2+2-1 = 3
- i=3: 2+3-1 = 4
Positions of 1s: {1, 3, 4}
Binary: 1 0 1 1 0

## μ = (1,1) — left edge (read bottom to top)
μ = (1,1,0) padded to length n=3:
- i=1: 2+1-1 = 2
- i=2: 2+2-1 = 3
- i=3: 2+3-0 = 5
Positions of 1s: {2, 3, 5}
Binary: 0 1 1 0 1

## ν = (1,1) — right edge (read top to bottom)
Same as μ: 0 1 1 0 1

## Triangle layout

Molev's convention: λ on bottom (left→right), λ goes up-left, μ goes up-right.
Wait — from the Molev slide, the triangle has:
- λ on the left edge (bottom-left to top, reading upward)
- μ on the right edge (top to bottom-right, reading downward)
- ν on the bottom edge (left to right)

Actually, looking at Molev's slide again:
- The triangle has arrows: λ↗ on left, μ↘ on right, ν→ on bottom
- Bottom edge: ν = 1010 (left to right)
- Left edge: λ = 1001 (bottom to top)
- Right edge: μ = 0101 (top to bottom)

So the convention is c^ν_{λμ}, not c^λ_{μν}! Molev uses c^ν_{λμ}.

For our problem: c^λ_{μν} with λ=(2,1,1), μ=ν=(1,1).

Using Molev's triangle layout (c^ν_{λμ}), we'd put:
- Left edge (↗): λ = first argument = μ = (1,1) → 01101
- Right edge (↘): μ = second argument = ν = (1,1) → 01101
- Bottom edge (→): ν = result = λ = (2,1,1) → 10110

## Binary strings on edges of size-5 triangle

Bottom (→): 1 0 1 1 0
Left (↗):   0 1 1 0 1
Right (↘):  0 1 1 0 1

## Internal structure

The size-5 triangle has 25 unit triangles (upward) and some downward triangles,
plus internal edges to label.

Rows from bottom to top:
- Row 1 (bottom): 5 up-triangles, 4 down-triangles, 9 internal horizontal edges
- Row 2: 3 up-triangles, 2 down-triangles
- Row 3: 3 up-triangles, 2 down-triangles

Actually for side length N=5:
- Row 1 (bottom): 9 unit triangles (5 up + 4 down)
- Row 2: 7 unit triangles (4 up + 3 down)
- Row 3: 5 (3 up + 2 down)
- Row 4: 3 (2 up + 1 down)
- Row 5 (top): 1 (1 up + 0 down)
Total = 9+7+5+3+1 = 25 triangles (15 up + 10 down)

This is getting complex. Let me try the ZJ example first, which has a SMALLER triangle.
Wait — what size triangle did ZJ use?

ZJ says n = N+(λ) - N-(λ). For λ=(2,1,1):
In the Fock space picture, the state |λ⟩ has particles at positions...
N+(λ) = rightmost position of a + edge, N-(λ) = leftmost position of a - edge.
For λ=(2,1,1) encoded as ...---++--+++...
The - positions (from the infinite encoding) are at: ..., -3, -2, 0, ...
Wait, this is getting into ZJ's infinite Fock space formalism.

But ZJ says k = height of λ = 3, and n = N+ - N- for the triangle.
For the finite version, n = side length = n_rows + n_cols of bounding rectangle.

Actually, I think ZJ and Molev use the SAME triangle — it's N=5 for our example.

This is a big puzzle (25 tiles). Let me try solving it computationally instead of by hand.
