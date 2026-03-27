# Binary Encoding — Working Notes

## Calibrating against Molev's example

Molev slide: λ=(2) inside 2×2 rectangle → binary code `1001`

The 2×2 rectangle with λ=(2):
```
+--+--+
|  |  |   row 2 (top): empty
+--+--+
|##|##|   row 1 (bottom): 2 boxes
+--+--+
```

### Trying: boundary path of the COMPLEMENT

The complement of λ=(2) inside 2×2 is the partition obtained by
rotating the empty boxes 180°. Empty boxes are the entire top row = (2).
So complement = (2).

Complement boundary from bottom-left, 1=up 0=right:
Path traces (2) = right,right,up,up = 0,0,1,1 → `0011`. Not `1001`.

### Trying: encode via column heights

λ=(2) in 2×2: column heights are (1,1).
For each column position j=1..m, put a 1 at position (n-λ'_j) where λ' is conjugate.
λ'=(1,1). Positions where there's a "step down" in the boundary...

### Trying: NE boundary, bottom to top, 1=right 0=up

λ=(2) in 2×2. The NE boundary of λ goes:
Start at bottom-right of row 1 (since row 1 has 2 boxes, we're at (2,0)).
Go up to (2,1). Then left to (0,1)... no, the NE boundary is the RIGHT+TOP edge.

### Trying: the KTW encoding directly

From Knutson-Tao-Woodward: the binary string for partition λ inside n×m
has length n+m. Place 1s at positions {m+1-λ_1, m+2-λ_2, ..., m+n-λ_n}
(using 1-indexed positions from left).

Wait no. Standard KTW: the string has n 1s and m 0s. The 1s are at positions
{λ_i + n+1-i : i=1..n} ... let me try the simpler version.

For λ inside n×m, binary string b of length n+m:
- b encodes a path from (0,0) to (m,n) where 0=right, 1=up
- The partition is the region between this path and the "staircase" path

Actually the simplest: the string has 0s and 1s, with n 1s and m 0s.
Position the 1s at indices {λ_i + n - i + 1} for i=1..n, using 1-based indexing...
no wait, I think position of 1s = {m - λ_i + i} for i = 1..n.

### Let me just try all possibilities for λ=(2) in 2×2

n=2, m=2, length=4, two 1s and two 0s.
Possible strings: 0011, 0101, 0110, 1001, 1010, 1100.

Molev says λ=(2) → 1001, μ=(1) → 0101, ν=(2,1) → 1010.

For λ=(2): 1s at positions 1,4 (1-indexed).
For μ=(1): 1s at positions 2,4.
For ν=(2,1): 1s at positions 1,3.

Pattern for λ_i in n=2, m=2: positions of 1s = {m+1-λ_i+i-1} = {3-λ_i+i-1} = {2-λ_i+i}?
- λ=(2): i=1: 2-2+1=1 ✓, i=2: 2-0+2=4 ✓  → positions {1,4} → 1001 ✓✓✓
- μ=(1): i=1: 2-1+1=2, i=2: 2-0+2=4 → positions {2,4} → 0101 ✓✓✓
- ν=(2,1): i=1: 2-2+1=1, i=2: 2-1+2=3 → positions {1,3} → 1010 ✓✓✓

## FORMULA: position of i-th 1 = m + i - λ_i  (1-indexed, i=1..n)

Or equivalently: b_j = 1 if j ∈ {m+i-λ_i : i=1..n}, else b_j = 0.
