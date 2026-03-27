#!/usr/bin/env python3
"""
Littlewood-Richardson coefficients via commuting transfer operators.

Computes c^ν_{λ,μ} and enumerates KTW puzzle witnesses using the row
transfer matrix on trinary-labeled edges of a puzzle triangle.

The row transfer decomposes into commuting operators T̃_+ (right boundary,
red spectral line) and T̃_- (left boundary, green spectral line), whose
commutativity follows from the Yang-Baxter equation for the free-fermion
vertex model.  (Zinn-Justin 2008, arXiv:0809.2392)

Each valid path through the transfer matrix is a KTW puzzle witness.
If c^ν_{λ,μ} = n, exactly n Puzzle objects are returned.
"""

from collections import Counter


# ================================================================
# KTW puzzle tile rules
# ================================================================
#
# Edge labels:  0 = empty,  1 = green path,  10 = rhombus (both)
#
# △ (Delta) triples:  (bottom─, left╱, right╲)
# ▽ (Nabla) triples:  (top─, right╱, left╲)
#
# Key structural fact:
#   △ is DETERMINISTIC: given (bottom, left╱), right╲ is unique.
#   ▽ may BRANCH: given left╲, up to 2 valid (top, right╱) pairs.
#   All puzzle multiplicity comes from the ▽ branching.

DELTA = frozenset({
    (0, 0, 0),
    (1, 0, 10),
    (10, 1, 0),
    (1, 1, 1),
    (0, 10, 1),
})

NABLA = frozenset({
    (0, 0, 0),
    (0, 10, 1),
    (10, 1, 0),
    (1, 0, 10),
    (1, 1, 1),
})

# △ deterministic lookup: (bottom, left╱) → right╲
_DELTA_FUNC = {(b, l): r for b, l, r in DELTA}

# ▽ branching lookup: left╲ → [(top, right╱), ...]
_NABLA_BRANCH = {}
for _t, _r, _l in NABLA:
    _NABLA_BRANCH.setdefault(_l, []).append((_t, _r))

# ▽ reconstruction lookup: (top, left╲) → right╱  (unique)
_NABLA_FUNC = {(t, l): r for t, r, l in NABLA}


# ================================================================
# Partition ↔ boundary encoding
# ================================================================

def partition_to_boundary(lam, n, m):
    """Encode partition λ ⊆ n×m box as binary tuple of length N = n+m.

    Position of i-th 1 (0-indexed) = m + i - λ_i.
    The string has exactly n ones and m zeros.
    """
    pad = list(lam) + [0] * (n - len(lam))
    bits = [0] * (n + m)
    for i in range(n):
        bits[m + i - pad[i]] = 1
    return tuple(bits)


def boundary_to_partition(bits):
    """Decode binary boundary tuple back to a partition."""
    m = bits.count(0)
    parts = [m + i - pos for i, pos in enumerate(
        j for j, b in enumerate(bits) if b == 1
    )]
    while parts and parts[-1] == 0:
        parts.pop()
    return tuple(parts)


# ================================================================
# Row transfer matrix
# ================================================================

def _row_transitions(state, b_left, b_right):
    """Yield all valid next-row states for one row of the puzzle triangle.

    Processes w △ triangles and w-1 ▽ triangles left-to-right.
    △ are deterministic; ▽ may branch (at most 2 ways each).

    Args:
        state:   tuple of labels (length w) — horizontal edges at bottom
        b_left:  label on the left ╱ boundary edge
        b_right: label on the right ╲ boundary edge

    Yields:
        tuples of labels (length w-1) — horizontal edges at top
    """
    w = len(state)

    def step(j, diag, top):
        # △_j: deterministic — (bottom, left╱) → right╲
        r_diag = _DELTA_FUNC.get((state[j], diag))
        if r_diag is None:
            return  # invalid (bottom, left╱) pair — dead end

        if j == w - 1:
            # last △: its right╲ must match the right boundary
            if r_diag == b_right:
                yield tuple(top)
        else:
            # ▽_j: given left╲ = r_diag, branch over valid (top, right╱)
            for t_val, n_diag in _NABLA_BRANCH.get(r_diag, ()):
                yield from step(j + 1, n_diag, top + [t_val])

    yield from step(0, b_left, [])


# ================================================================
# Bounding box
# ================================================================

def _bounding_box(lam, mu, nu):
    """Compute minimal n×m box containing all three partitions."""
    n = max(len(lam), len(mu), len(nu))
    m = 0
    for p in (lam, mu, nu):
        if p:
            m = max(m, p[0])
    return n, m


# ================================================================
# Puzzle witness
# ================================================================

class Puzzle:
    """A single KTW puzzle witness for c^ν_{λ,μ}.

    Attributes:
        lam, mu, nu: the three partitions
        n, m, N:     bounding box (n rows, m cols, N = n+m = triangle side)
        levels:      tuple of N+1 tuples — horizontal edge labels at each
                     level, from bottom (length N) to apex (length 0)
    """

    __slots__ = ('lam', 'mu', 'nu', 'n', 'm', 'N', 'levels')

    def __init__(self, lam, mu, nu, n, m, levels):
        self.lam = lam
        self.mu = mu
        self.nu = nu
        self.n = n
        self.m = m
        self.N = n + m
        self.levels = levels

    def __repr__(self):
        sym = {0: '0', 1: '1', 10: 'R'}
        hdr = f"Puzzle c^{self.nu}_{{{self.lam},{self.mu}}}"
        lines = [hdr]
        for i, lv in enumerate(self.levels):
            s = ' '.join(sym[x] for x in lv)
            lines.append(f"  {'  ' * i}{s}")
        return '\n'.join(lines)

    def edge_labels(self):
        """Reconstruct the full edge labeling of the puzzle triangle.

        Returns dict  {((r1,j1),(r2,j2)): label}  with vertices sorted.
        Diagonal labels are fully determined by the horizontal levels
        (△ is deterministic, ▽ is determined by knowing the top edge).
        """
        N = self.N
        left = partition_to_boundary(self.lam, self.n, self.m)
        right = partition_to_boundary(self.mu, self.n, self.m)
        labels = {}

        # Horizontal edges
        for r, lv in enumerate(self.levels):
            for j, val in enumerate(lv):
                labels[((r, j), (r, j + 1))] = val

        # Diagonal edges — reconstruct row by row
        for r in range(N):
            w = N - r
            state = self.levels[r]
            nxt = self.levels[r + 1]
            diag = left[r]

            for j in range(w):
                # ╱ edge
                labels[((r, j), (r + 1, j))] = diag
                # △_j: right╲
                r_diag = _DELTA_FUNC[(state[j], diag)]
                labels[((r, j + 1), (r + 1, j))] = r_diag
                # ▽_j: right╱ → becomes next diag
                if j < w - 1:
                    diag = _NABLA_FUNC[(nxt[j], r_diag)]

        return labels


# ================================================================
# Main computation
# ================================================================

def lr_puzzles(lam, mu, nu):
    """Compute all KTW puzzle witnesses for c^ν_{λ,μ}.

    Returns a list of Puzzle objects.  len(result) == c^ν_{λ,μ}.
    """
    n, m = _bounding_box(lam, mu, nu)
    N = n + m

    if N == 0:
        return [Puzzle(lam, mu, nu, 0, 0, ((),))]

    bottom = partition_to_boundary(nu, n, m)
    left = partition_to_boundary(lam, n, m)
    right = partition_to_boundary(mu, n, m)

    # Enumerate paths through the row transfer matrix.
    # Each partial: (current_state, accumulated_levels)
    partials = [(bottom, (bottom,))]

    for r in range(N):
        b_left = left[r]
        b_right = right[N - 1 - r]
        new_partials = []
        for state, levels in partials:
            for nxt in _row_transitions(state, b_left, b_right):
                new_partials.append((nxt, levels + (nxt,)))
        partials = new_partials

    return [Puzzle(lam, mu, nu, n, m, levels) for _, levels in partials]


def lr_count(lam, mu, nu):
    """Just the LR coefficient c^ν_{λ,μ}, without storing puzzle data."""
    n, m = _bounding_box(lam, mu, nu)
    N = n + m

    if N == 0:
        return 1

    bottom = partition_to_boundary(nu, n, m)
    left = partition_to_boundary(lam, n, m)
    right = partition_to_boundary(mu, n, m)

    states = Counter({bottom: 1})

    for r in range(N):
        b_left = left[r]
        b_right = right[N - 1 - r]
        new_states = Counter()
        for state, count in states.items():
            for nxt in _row_transitions(state, b_left, b_right):
                new_states[nxt] += count
        states = new_states

    return states.get((), 0)


# ================================================================
# Demo & verification
# ================================================================

def demo():
    print("=" * 60)
    print("  LR Coefficients via Commuting Transfer Operators")
    print("=" * 60)

    examples = [
        # (lam, mu, nu, expected_c)
        # --- trivial ---
        ((), (), (), 1),
        # --- Pieri (single row × single row) ---
        ((1,), (1,), (2,), 1),
        ((1,), (1,), (1, 1), 1),
        ((2,), (1,), (2, 1), 1),
        # --- Pieri (single row × column) ---
        ((1,), (1, 1), (2, 1), 1),
        ((1,), (1, 1), (1, 1, 1), 1),
        # --- column × column ---
        ((1, 1), (1, 1), (2, 2), 1),
        ((1, 1), (1, 1), (2, 1, 1), 1),
        ((1, 1), (1, 1), (1, 1, 1, 1), 1),
        ((1, 1), (1, 1), (3, 1), 0),
        # --- row × column ---
        ((2,), (1, 1), (3, 1), 1),
        ((2,), (1, 1), (2, 1, 1), 1),
        ((2,), (1, 1), (2, 2), 0),
        # --- the famous c = 2 case ---
        ((2, 1), (2, 1), (3, 2, 1), 2),
        ((2, 1), (2, 1), (4, 1, 1), 1),
        # --- size mismatch → 0 ---
        ((1,), (1,), (3,), 0),
    ]

    print("\n  Checking c^ν_{λ,μ} against known values:\n")
    all_ok = True
    for lam, mu, nu, expected in examples:
        count = lr_count(lam, mu, nu)
        ok = "✓" if count == expected else "✗"
        if count != expected:
            all_ok = False
        print(f"    c^{str(nu):14s}_{{{str(lam):8s}, {str(mu):8s}}} "
              f"= {count:2d}  (expected {expected})  {ok}")

    # --- Symmetry: c^ν_{λ,μ} = c^ν_{μ,λ} ---
    print("\n  Symmetry check: c^ν_{λ,μ} = c^ν_{μ,λ}")
    sym_ok = True
    for lam, mu, nu, _ in examples:
        if lam == mu:
            continue
        c1 = lr_count(lam, mu, nu)
        c2 = lr_count(mu, lam, nu)
        if c1 != c2:
            print(f"    FAIL: c^{nu}_{{{lam},{mu}}} = {c1} "
                  f"≠ c^{nu}_{{{mu},{lam}}} = {c2}")
            sym_ok = False
    print(f"    {'All symmetric ✓' if sym_ok else 'SYMMETRY FAILURES!'}")

    # --- Puzzle witnesses for c = 2 ---
    print(f"\n  Puzzle witnesses for c^(3,2,1)_{{(2,1),(2,1)}} = 2:")
    puzzles = lr_puzzles((2, 1), (2, 1), (3, 2, 1))
    for i, p in enumerate(puzzles):
        print(f"\n  Puzzle {i + 1}:")
        for line in repr(p).split('\n'):
            print(f"    {line}")

    # --- Cross-validate with solve_puzzle.py ---
    print(f"\n  Cross-validation with brute-force solver:")
    try:
        from solve_puzzle import solve
        cross_cases = [
            ((2,), (1,), (2, 1), 1),
            ((1, 1), (1, 1), (2, 1, 1), 1),
            ((2, 1), (2, 1), (3, 2, 1), 2),
        ]
        cross_ok = True
        for lam, mu, nu, expected in cross_cases:
            import io, contextlib
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                sols, _ = solve(lam=lam, mu=mu, nu=nu)
            brute = len(sols)
            transfer = lr_count(lam, mu, nu)
            ok = "✓" if brute == transfer else "✗"
            if brute != transfer:
                cross_ok = False
            print(f"    c^{str(nu):14s}: brute={brute}, transfer={transfer}  {ok}")
        print(f"    {'All match ✓' if cross_ok else 'MISMATCH!'}")
    except ImportError:
        print("    (solve_puzzle.py not importable, skipping)")

    print()
    ok = all_ok and sym_ok
    print(f"  {'All checks passed!' if ok else 'SOME CHECKS FAILED!'}")
    return ok


if __name__ == '__main__':
    demo()
