#!/usr/bin/env python3
"""
Transfer matrix operators from Zinn-Justin (2008)
"Littlewood-Richardson coefficients and integrable tilings", arXiv:0809.2392

Implements operators on partitions represented as binary strings
in the fermionic Fock space F and the tiling Fock space G.

ENCODING CONVENTION (Section 3.1 + 5.2 of the paper):
  Partition λ inside k × m box  →  binary string of length n = k + m
  Trace the NE boundary of the Young diagram from (0,0) to (m,k):
    1 ≡ '−' = vertical step (↑) = green particle
    0 ≡ '+' = horizontal step (→) = red particle

  Example: λ = (3,2,1) in 3×4 box
    Boundary: →↑→↑→↑→ = 0101010
    Binary:   (0,1,0,1,0,1,0)

OPERATORS (Sections 3.4, 3.6):
  T_free(u)|λ⟩  = Σ_{μ⊆λ, λ/μ h-strip} u^{|λ/μ|} |μ⟩   removes horizontal strips
  T̄_free(u)|λ⟩ = Σ_{μ⊆λ, λ/μ v-strip} u^{|λ/μ|} |μ⟩   removes vertical strips
  T̃_+(u) on G_free: acts as T_free on ν  (right/red part)
  T̃_-(u) on G_free: acts as T̄_free on μ (left/green part)

KEY FORMULAS (Lemma 1, Theorem 4):
  s_{λ/μ}(u₁,…,uₙ) = ⟨μ| ∏ T_free(uᵢ) |λ⟩
  s_{λᵀ/μᵀ}(u₁,…,uₙ) = ⟨μ| ∏ T̄_free(uᵢ) |λ⟩
  s_λ(u/v) = Σ_{μ,ν} c^λ_{μ,ν} s_ν(u) s_{μᵀ}(v)     (coproduct = LR coefficients)
"""

from collections import defaultdict


# ================================================================
# Part I: Partition ↔ Binary string
# ================================================================

def partition_to_binary(lam, k, m):
    """
    Encode partition λ in a k×m box as a binary string of length k+m.

    Traces the NE boundary from (0,0) to (m,k):
      1 = vertical step (↑) = green particle
      0 = horizontal step (→) = red particle

    >>> partition_to_binary((2,1), 2, 3)
    (0, 1, 0, 1, 0)
    >>> partition_to_binary((), 2, 3)
    (1, 1, 0, 0, 0)
    >>> partition_to_binary((3,3), 2, 3)
    (0, 0, 0, 1, 1)
    """
    lam = list(lam) + [0] * (k - len(lam))
    assert len(lam) == k
    assert all(0 <= p <= m for p in lam), f"Parts must be in [0,{m}], got {lam}"
    assert all(lam[i] >= lam[i + 1] for i in range(k - 1)), f"Not weakly decreasing: {lam}"

    bits = []
    prev_width = 0
    for i in range(k - 1, -1, -1):  # bottom row (shortest) to top row (longest)
        bits.extend([0] * (lam[i] - prev_width))  # horizontal steps
        bits.append(1)  # vertical step
        prev_width = lam[i]
    bits.extend([0] * (m - prev_width))  # remaining horizontal steps
    return tuple(bits)


def binary_to_partition(b):
    """
    Decode binary string to partition. Infers k and m from the string.

    >>> binary_to_partition((0, 1, 0, 1, 0))
    (2, 1)
    >>> binary_to_partition((1, 1, 0, 0, 0))
    (0, 0)
    """
    parts = []
    col = 0
    for bit in b:
        if bit == 0:
            col += 1
        else:
            parts.append(col)
    parts.reverse()
    # Strip trailing zeros for cleaner output
    while parts and parts[-1] == 0:
        parts.pop()
    return tuple(parts) if parts else ()


def fmt_binary(b):
    """Pretty-print binary string using ±/green-red notation."""
    return ''.join('−' if x == 1 else '+' for x in b)


def partition_size(lam):
    return sum(lam)


def conjugate(lam):
    """Conjugate (transpose) of a partition."""
    if not lam:
        return ()
    return tuple(sum(1 for p in lam if p > j) for j in range(lam[0]))


# ================================================================
# Part II: Strip conditions
# ================================================================

def is_horizontal_strip(lam, mu):
    """
    Check if λ/μ is a horizontal strip (at most one box per column).

    Equivalent to the INTERLACING condition:
      λ₁ ≥ μ₁ ≥ λ₂ ≥ μ₂ ≥ λ₃ ≥ … ≥ 0

    >>> is_horizontal_strip((3,2,1), (2,2,0))
    True
    >>> is_horizontal_strip((3,2,1), (3,1,1))
    False
    """
    k = max(len(lam), len(mu))
    L = list(lam) + [0] * (k - len(lam))
    M = list(mu) + [0] * (k - len(mu))

    for i in range(k):
        lower = L[i + 1] if i + 1 < k else 0
        if not (L[i] >= M[i] >= lower):
            return False
    return True


def is_vertical_strip(lam, mu):
    """
    Check if λ/μ is a vertical strip (at most one box per row).

    Equivalent to: μᵢ ∈ {λᵢ, λᵢ − 1} for all i, and μ weakly decreasing.

    >>> is_vertical_strip((3,2,1), (2,2,1))
    True
    >>> is_vertical_strip((3,2,1), (1,2,1))
    False
    """
    k = max(len(lam), len(mu))
    L = list(lam) + [0] * (k - len(lam))
    M = list(mu) + [0] * (k - len(mu))

    for i in range(k):
        if L[i] - M[i] not in (0, 1):
            return False
    return all(M[i] >= M[i + 1] for i in range(k - 1))


# ================================================================
# Part III: T_free(u) — removes horizontal strips
# ================================================================

def T_free(lam, k, m):
    """
    Compute T_free(u)|λ⟩ = Σ_μ u^{|λ/μ|} |μ⟩
    where the sum is over all μ ⊆ λ such that λ/μ is a horizontal strip.

    Binary string picture: red particles (+/0) cascade RIGHT past green (−/1).
    Each green-red crossing contributes one power of u.

    Returns: dict {μ_tuple: exponent_of_u}

    >>> sorted(T_free((2,1), 2, 3).items())
    [((1,), 2), ((1, 1), 1), ((2,), 1), ((2, 1), 0)]
    """
    L = list(lam) + [0] * (k - len(lam))
    results = {}

    def recurse(i, upper, mu_so_far):
        if i == k:
            mu = tuple(mu_so_far)
            # Strip trailing zeros
            while mu and mu[-1] == 0:
                mu = mu[:-1]
            power = partition_size(L) - partition_size(mu)
            results[mu] = power
            return
        # Interlacing: L[i] ≥ μ_i ≥ L[i+1], and μ_i ≤ upper (weakly decreasing)
        lower = L[i + 1] if i + 1 < k else 0
        top = min(L[i], upper)
        for mu_i in range(top, lower - 1, -1):
            recurse(i + 1, mu_i, mu_so_far + [mu_i])

    recurse(0, m, [])
    return results


def T_bar_free(lam, k, m):
    """
    Compute T̄_free(u)|λ⟩ = Σ_μ u^{|λ/μ|} |μ⟩
    where the sum is over all μ ⊆ λ such that λ/μ is a vertical strip.

    Binary string picture: green particles (−/1) each move LEFT by 0 or 1 step
    past a red (+/0). Each such move contributes one power of u.

    Returns: dict {μ_tuple: exponent_of_u}

    >>> sorted(T_bar_free((2,1), 2, 3).items())
    [((1,), 2), ((1, 1), 1), ((2,), 1), ((2, 1), 0)]
    """
    L = list(lam) + [0] * (k - len(lam))
    results = {}

    def recurse(i, upper, mu_so_far):
        if i == k:
            mu = tuple(mu_so_far)
            while mu and mu[-1] == 0:
                mu = mu[:-1]
            power = partition_size(L) - partition_size(mu)
            results[mu] = power
            return
        # Vertical strip: μ_i ∈ {L[i], L[i]-1}, also μ_i ≤ upper and μ_i ≥ 0
        for delta in (0, 1):
            mu_i = L[i] - delta
            if mu_i < 0 or mu_i > upper:
                continue
            recurse(i + 1, mu_i, mu_so_far + [mu_i])

    recurse(0, m, [])
    return results


# ================================================================
# Part IV: Operator algebra — polynomials in u, composition
# ================================================================

# A "state" is dict { partition_tuple : polynomial }
# A polynomial is dict { power_of_u : integer_coefficient }

def state_single(lam):
    """State |λ⟩ with coefficient 1."""
    return {lam: {0: 1}}


def apply_op(op_func, state, k, m):
    """
    Apply operator (T_free or T_bar_free) to a state.

    state: dict {partition: {power: coeff}}
    Returns: new state dict.
    """
    new_state = defaultdict(lambda: defaultdict(int))
    for lam, poly in state.items():
        targets = op_func(lam, k, m)
        for mu, extra_power in targets.items():
            for p, c in poly.items():
                new_state[mu][p + extra_power] += c

    # Clean zeros
    result = {}
    for mu, poly in new_state.items():
        cleaned = {p: c for p, c in poly.items() if c != 0}
        if cleaned:
            result[mu] = cleaned
    return result


def inner_product(bra, state):
    """⟨bra|state⟩ → polynomial {power: coeff}."""
    return dict(state.get(bra, {}))


def fmt_poly(poly, var='u'):
    """Pretty-print polynomial."""
    if not poly:
        return '0'
    terms = []
    for p in sorted(poly.keys()):
        c = poly[p]
        if c == 0:
            continue
        if p == 0:
            terms.append(str(c))
        elif c == 1:
            terms.append(var if p == 1 else f'{var}^{p}')
        else:
            terms.append(f'{c}{var}' if p == 1 else f'{c}{var}^{p}')
    return ' + '.join(terms) if terms else '0'


# ================================================================
# Part V: T̃_+ and T̃_- on G_free
# ================================================================

def T_tilde_plus(mu, nu, k1, m1, k2, m2):
    """
    T̃_+(u) on G_free: acts as T_free(u) on ν (right part), μ unchanged.

    |μ⟩ ⊔ |ν⟩  →  Σ u^power |μ⟩ ⊔ |ν'⟩

    Returns: dict { (μ, ν') : power_of_u }
    """
    targets = T_free(nu, k2, m2)
    return {(mu, nu_p): power for nu_p, power in targets.items()}


def T_tilde_minus(mu, nu, k1, m1, k2, m2):
    """
    T̃_-(u) on G_free: acts as T̄_free(u) on μ (left part), ν unchanged.

    |μ⟩ ⊔ |ν⟩  →  Σ u^power |μ'⟩ ⊔ |ν⟩

    Returns: dict { (μ', ν) : power_of_u }
    """
    targets = T_bar_free(mu, k1, m1)
    return {(mu_p, nu): power for mu_p, power in targets.items()}


# ================================================================
# Part VI: Verification
# ================================================================

def schur_via_transfer(lam, n_vars, k, m):
    """
    Compute s_λ(u,...,u) [n_vars copies] via ⟨∅| ∏ T_free(u) |λ⟩.
    Returns polynomial in u.
    """
    empty = ()
    state = state_single(lam)
    for _ in range(n_vars):
        state = apply_op(T_free, state, k, m)
    return inner_product(empty, state)


def schur_transpose_via_transfer(lam, n_vars, k, m):
    """
    Compute s_{λᵀ}(u,...,u) via ⟨∅| ∏ T̄_free(u) |λ⟩.
    """
    empty = ()
    state = state_single(lam)
    for _ in range(n_vars):
        state = apply_op(T_bar_free, state, k, m)
    return inner_product(empty, state)


def hook_content_formula(lam, n):
    """
    s_λ(1,...,1) [n copies] via hook-content formula.
    = ∏_{(i,j)∈λ} (n + j - i) / h(i,j)

    Must accumulate numerator and denominator separately to avoid
    premature integer division.
    """
    if not lam:
        return 1
    k = len(lam)
    lam_conj = conjugate(lam)
    numer = 1
    denom = 1
    for i in range(k):
        for j in range(lam[i]):
            content = j - i
            arm = lam[i] - j - 1
            leg = lam_conj[j] - i - 1
            hook = arm + leg + 1
            numer *= (n + content)
            denom *= hook
    assert numer % denom == 0, f"Not integer: {numer}/{denom}"
    return numer // denom


# ================================================================
# Part VII: Coproduct — Δ(s_λ) via transfer operators
# ================================================================
#
# The Hopf algebra coproduct on symmetric functions:
#   Δ(s_λ) = Σ_{μ,ν} c^λ_{μ,ν} s_μ ⊗ s_ν
#
# Method: decompose s_{λ/μ} = Σ_ν c^λ_{μ,ν} s_ν using the Kostka
# matrix.  Compute skew Kostka numbers K̃_{λ/μ,ρ} (= #SSYT of shape
# λ/μ and content ρ) via T_free chains, then invert the unitriangular
# Kostka matrix to extract individual c^λ_{μ,ν}.
#
# Why not evaluate at x=(1,...,1)?  The hook-content polynomials
# dim_n(ν) are linearly DEPENDENT for some partition sets — e.g.,
# dim_n(3,1) − 3·dim_n(2,2) + dim_n(2,1,1) = 0 identically in n.
# So the principal specialization cannot separate all Schur functions.


def _partitions_of(d, max_len, max_part):
    """Enumerate all partitions of d with ≤ max_len parts, each ≤ max_part."""
    if d == 0:
        yield ()
        return
    if max_len == 0 or max_part == 0:
        return
    for first in range(min(d, max_part), 0, -1):
        for rest in _partitions_of(d - first, max_len - 1, first):
            yield (first,) + rest


def _chain_count(start, strip_sizes, k, m):
    """Count endpoints of horizontal-strip chains from `start`.

    Removes strips of sizes strip_sizes[-1], strip_sizes[-2], ..., strip_sizes[0]
    (last part first, matching SSYT convention where highest entry fills last).

    Returns: dict {endpoint_partition: count}.
    """
    current = {start: 1}
    for sz in reversed(strip_sizes):
        nxt = {}
        for eta, cnt in current.items():
            for tau, power in T_free(eta, k, m).items():
                if power == sz:
                    nxt[tau] = nxt.get(tau, 0) + cnt
        current = nxt
    return current


def lr_via_transfer(lam, k_total=None, m_total=None):
    """
    Compute the coproduct Δ(s_λ) = Σ_{μ,ν} c^λ_{μ,ν} s_μ ⊗ s_ν
    via transfer matrix operators.

    Decomposes s_{λ/μ} into Schur functions by computing skew Kostka
    numbers K̃_{λ/μ,ρ} (via T_free chains) and inverting the Kostka
    matrix (unitriangular in dominance order) to extract c^λ_{μ,ν}.

    Returns: dict {(μ, ν): c^λ_{μ,ν}} for all pairs with c > 0.
    """
    lam = tuple(p for p in lam if p > 0) if lam else ()
    total = partition_size(lam)
    if total == 0:
        return {((), ()): 1}

    k = k_total or len(lam)
    m = m_total or lam[0]

    # Candidate ν grouped by size d = |ν|, sorted by decreasing dominance
    candidate_groups = {}
    for d in range(total + 1):
        cands = list(_partitions_of(d, k, m))
        if cands:
            # Sort key: negative partial sums (most dominant first)
            ml = max(len(c) for c in cands)
            def _key(p, _ml=ml):
                s, out = 0, []
                for x in list(p) + [0] * (_ml - len(p)):
                    s += x
                    out.append(-s)
                return out
            candidate_groups[d] = sorted(cands, key=_key)

    # Precompute Kostka matrices K[d][i][j] = K_{cands_d[i], cands_d[j]}
    # K is unitriangular: K[i][i] = 1, K[i][j] = 0 when cands[i] < cands[j]
    kostka = {}
    for d, cands in candidate_groups.items():
        nc = len(cands)
        mat = [[0] * nc for _ in range(nc)]
        for i in range(nc):
            for j in range(nc):
                endpoints = _chain_count(cands[i], cands[j], k, m)
                mat[i][j] = endpoints.get((), 0)
        kostka[d] = mat

    # Compute skew Kostka: for each ρ, {μ: K̃_{λ/μ,ρ}} for all μ at once
    skew_kostka = {}   # (d, j) → {μ: K̃ value}
    all_mu = set()
    for d, cands in candidate_groups.items():
        for j, rho in enumerate(cands):
            endpoints = _chain_count(lam, rho, k, m)
            skew_kostka[(d, j)] = endpoints
            all_mu.update(endpoints.keys())

    # For each μ, invert Kostka matrix to get c^λ_{μ,ν}
    result = {}
    for mu in all_mu:
        d = total - partition_size(mu)
        if d not in candidate_groups:
            continue
        cands = candidate_groups[d]
        nc = len(cands)
        K = kostka[d]

        # sk[i] = K̃_{λ/μ, cands[i]}
        sk = [skew_kostka.get((d, i), {}).get(mu, 0) for i in range(nc)]
        if all(v == 0 for v in sk):
            continue

        # Forward substitution on K^T c = sk  (K unitriangular)
        c = [0] * nc
        for i in range(nc):
            c[i] = sk[i] - sum(c[j] * K[j][i] for j in range(i))

        for i in range(nc):
            if c[i] > 0:
                result[(mu, cands[i])] = c[i]

    return result


# ================================================================
# Demo
# ================================================================

def demo():
    print("=" * 70)
    print("  Zinn-Justin Transfer Matrix Operators (arXiv:0809.2392)")
    print("=" * 70)

    # --- Binary string encoding ---
    k, m = 3, 4
    print(f"\n[1] Binary string encoding (k={k}, m={m} box)")
    print("-" * 50)
    examples = [(3, 2, 1), (4, 2, 0), (4, 4, 4), (0, 0, 0), (2, 1, 0)]
    for lam in examples:
        b = partition_to_binary(lam, k, m)
        lam_back = binary_to_partition(b)
        print(f"  λ = {str(lam):12s}  binary: {fmt_binary(b):10s}  = {b}")
        assert lam_back == tuple(p for p in lam if p > 0) or (lam_back == () and all(p == 0 for p in lam))

    # --- T_free: horizontal strip removal ---
    lam = (3, 2, 1)
    print(f"\n[2] T_free(u)|{lam}⟩  —  removes horizontal strips")
    print("-" * 50)
    results = T_free(lam, k, m)
    b_lam = partition_to_binary(lam, k, m)
    for mu in sorted(results, key=lambda x: (results[x], x)):
        power = results[mu]
        # Pad mu to k parts for binary conversion
        mu_padded = mu + (0,) * (k - len(mu))
        b_mu = partition_to_binary(mu_padded, k, m)
        print(f"  u^{power}  |{str(mu):12s}⟩   {fmt_binary(b_lam)} → {fmt_binary(b_mu)}"
              f"   λ/μ = {partition_size(lam) - partition_size(mu)} boxes")
    print(f"  Total: {len(results)} terms")

    # Verify horizontal strip condition
    for mu, _ in results.items():
        mu_padded = mu + (0,) * (k - len(mu))
        assert is_horizontal_strip(lam, mu_padded), f"FAIL: {lam}/{mu_padded} not h-strip"

    # --- T̄_free: vertical strip removal ---
    print(f"\n[3] T̄_free(u)|{lam}⟩  —  removes vertical strips")
    print("-" * 50)
    results_v = T_bar_free(lam, k, m)
    for mu in sorted(results_v, key=lambda x: (results_v[x], x)):
        power = results_v[mu]
        mu_padded = mu + (0,) * (k - len(mu))
        b_mu = partition_to_binary(mu_padded, k, m)
        print(f"  u^{power}  |{str(mu):12s}⟩   {fmt_binary(b_lam)} → {fmt_binary(b_mu)}"
              f"   λ/μ = {partition_size(lam) - partition_size(mu)} boxes")
    print(f"  Total: {len(results_v)} terms")

    for mu, _ in results_v.items():
        mu_padded = mu + (0,) * (k - len(mu))
        assert is_vertical_strip(lam, mu_padded), f"FAIL: {lam}/{mu_padded} not v-strip"

    # --- Difference between T_free and T̄_free ---
    print(f"\n[4] Comparing T_free vs T̄_free outputs")
    print("-" * 50)
    only_h = set(results) - set(results_v)
    only_v = set(results_v) - set(results)
    both = set(results) & set(results_v)
    print(f"  Only in T_free  (h-strip only): {len(only_h)}")
    for mu in sorted(only_h, key=lambda x: (results[x], x)):
        print(f"    μ = {mu}, removed {partition_size(lam) - partition_size(mu)} boxes")
    print(f"  Only in T̄_free (v-strip only): {len(only_v)}")
    for mu in sorted(only_v, key=lambda x: (results_v[x], x)):
        print(f"    μ = {mu}, removed {partition_size(lam) - partition_size(mu)} boxes")
    print(f"  In both (rook strips):          {len(both)}")

    # --- Schur function verification ---
    print(f"\n[5] Schur function verification: s_λ(u,...,u) = ⟨∅|∏ T_free(u)|λ⟩")
    print("-" * 50)
    test_cases = [(2, 1), (3, 1), (2, 2), (3, 2, 1)]
    for lam in test_cases:
        kk = len(lam)
        mm = lam[0]
        for n in [2, 3, 4]:
            poly = schur_via_transfer(lam, n, kk, mm)
            # s_λ(u,...,u) = u^{|λ|} * s_λ(1,...,1)
            # When n < height(λ), s_λ(1,...,1) = 0 if any hook-content factor is 0
            expected_count = hook_content_formula(lam, n)
            if expected_count == 0:
                expected_poly = {}
            else:
                expected_poly = {partition_size(lam): expected_count}
            ok = "✓" if poly == expected_poly else "✗"
            print(f"  s_{str(lam):10s}(u×{n}) = {fmt_poly(poly):15s}"
                  f"  expected {fmt_poly(expected_poly):15s}  {ok}")

    # --- T̄_free verification: s_{λ^T} ---
    print(f"\n[6] Transpose Schur: s_{{λᵀ}}(u,...,u) = ⟨∅|∏ T̄_free(u)|λ⟩")
    print("-" * 50)
    for lam in [(2, 1), (3, 1), (3, 2, 1)]:
        lam_T = conjugate(lam)
        kk = len(lam)
        mm = lam[0]
        for n in [2, 3]:
            poly = schur_transpose_via_transfer(lam, n, kk, mm)
            expected_count = hook_content_formula(lam_T, n)
            if expected_count == 0:
                expected_poly = {}
            else:
                expected_poly = {partition_size(lam): expected_count}
            ok = "✓" if poly == expected_poly else "✗"
            print(f"  s_{str(lam_T):10s}(u×{n}) via T̄|{lam}⟩ = {fmt_poly(poly):15s}"
                  f"  expected {fmt_poly(expected_poly):15s}  {ok}")

    # --- T̃ operators on G_free ---
    print(f"\n[7] T̃ operators on G_free = |μ⟩ ⊔ |ν⟩")
    print("-" * 50)
    mu, nu = (2, 1), (1,)
    k1, m1 = 2, 3
    k2, m2 = 1, 2
    print(f"  μ = {mu}  binary: {fmt_binary(partition_to_binary(mu, k1, m1))}")
    print(f"  ν = {nu}  binary: {fmt_binary(partition_to_binary(nu, k2, m2))}")

    print(f"\n  T̃_+(u): acts as T_free(u) on ν (right/red part)")
    tp = T_tilde_plus(mu, nu, k1, m1, k2, m2)
    for (m2_part, n2_part), power in sorted(tp.items(), key=lambda x: x[1]):
        print(f"    u^{power}  |{m2_part}⟩ ⊔ |{n2_part}⟩")

    print(f"\n  T̃_-(u): acts as T̄_free(u) on μ (left/green part)")
    tm = T_tilde_minus(mu, nu, k1, m1, k2, m2)
    for (m2_part, n2_part), power in sorted(tm.items(), key=lambda x: x[1]):
        print(f"    u^{power}  |{m2_part}⟩ ⊔ |{n2_part}⟩")

    # --- Binary string operation visualization ---
    print(f"\n[8] Binary string operations — detailed view")
    print("-" * 50)
    lam = (2, 1)
    k, m = 2, 3
    b_lam = partition_to_binary(lam, k, m)
    print(f"\n  λ = {lam}  in {k}×{m} box")
    print(f"  Binary: {fmt_binary(b_lam)} = {b_lam}")
    print(f"  Positions of '−' (green/1): {[i for i,x in enumerate(b_lam) if x==1]}")
    print(f"  Positions of '+' (red/0):   {[i for i,x in enumerate(b_lam) if x==0]}")

    print(f"\n  T_free: reds (+) cascade RIGHT past greens (−).")
    print(f"  Each crossing = one power of u. Result = horizontal strip removed.")
    results = T_free(lam, k, m)
    for mu in sorted(results, key=lambda x: (results[x], x)):
        power = results[mu]
        mu_padded = mu + (0,) * (k - len(mu))
        b_mu = partition_to_binary(mu_padded, k, m)
        # Show which positions changed
        changes = []
        for i in range(len(b_lam)):
            if b_lam[i] != b_mu[i]:
                changes.append(i)
        print(f"    {fmt_binary(b_lam)} → {fmt_binary(b_mu)}  u^{power}  "
              f"μ={mu}  changed at {changes}")

    print(f"\n  T̄_free: greens (−) move LEFT by ≤1 step past reds (+).")
    print(f"  Each move = one power of u. Result = vertical strip removed.")
    results_v = T_bar_free(lam, k, m)
    for mu in sorted(results_v, key=lambda x: (results_v[x], x)):
        power = results_v[mu]
        mu_padded = mu + (0,) * (k - len(mu))
        b_mu = partition_to_binary(mu_padded, k, m)
        changes = []
        for i in range(len(b_lam)):
            if b_lam[i] != b_mu[i]:
                changes.append(i)
        print(f"    {fmt_binary(b_lam)} → {fmt_binary(b_mu)}  u^{power}  "
              f"μ={mu}  changed at {changes}")

    # --- Coproduct Δ(s_λ) via transfer operators ---
    print(f"\n[9] Coproduct Δ(s_λ) = Σ c^λ_{{μ,ν}} s_μ ⊗ s_ν  (via transfer matrices)")
    print("-" * 50)
    for lam in [(2, 1), (3, 2, 1)]:
        kk, mm = len(lam), lam[0]
        coprod = lr_via_transfer(lam, kk, mm)
        print(f"\n  Δ(s_{lam}):  {len(coprod)} non-zero terms")
        for (mu, nu) in sorted(coprod):
            print(f"    c^{lam}_{{{mu}, {nu}}} = {coprod[(mu, nu)]}")
        try:
            from lr_coefficients import lr_count as _lr_count
            ok = all(_lr_count(mu, nu, lam) == c
                     for (mu, nu), c in coprod.items())
            print(f"    Puzzle cross-check: {'✓' if ok else 'FAIL'}")
        except ImportError:
            pass

    print()


if __name__ == '__main__':
    demo()
