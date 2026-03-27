"""
KTW Puzzle Solver — correct 3-label formulation.

Edge labels: 0, 1, 10 (where 10 is the rhombus label).
Boundary edges use only {0, 1}. Internal edges use {0, 1, 10}.

From SageMath's knutson_tao_puzzles.py:

Delta (△) triples (south─, NW╱, NE╲):
  (0, 0, 0)
  (0, 1, 10)
  (1, 10, 0)
  (1, 1, 1)
  (10, 0, 1)

Nabla (▽) triples (north─, SE╲, SW╱):
  (0, 0, 0)
  (0, 10, 1)
  (10, 1, 0)
  (1, 0, 10)
  (1, 1, 1)
"""

import math

def edge(v1, v2):
    return (min(v1, v2), max(v1, v2))

def partition_to_binary(rho, n, m):
    """Partition → binary string of length n+m.
    Position of i-th 1 (1-indexed) = m + i - rho[i-1]."""
    rho_padded = list(rho) + [0] * (n - len(rho))
    ones = set()
    for i in range(n):
        pos = m + (i + 1) - rho_padded[i]
        ones.add(pos - 1)
    return [1 if j in ones else 0 for j in range(n + m)]

# Allowed triples from SageMath.
# △ (Delta): (south─, NW╱, NE╲)
# My coordinates: △(r,j) has bottom─, left╱, right╲
# So: (bottom, left, right)
DELTA_TRIPLES = {
    (0, 0, 0),    # 0/0\0:  south=0, NW=0, NE=0
    (1, 0, 10),   # 0/1\10: south=1, NW=0, NE=10
    (10, 1, 0),   # 1/10\0: south=10, NW=1, NE=0
    (1, 1, 1),    # 1/1\1:  south=1, NW=1, NE=1
    (0, 10, 1),   # 10/0\1: south=0, NW=10, NE=1
}

# ▽ (Nabla): (north─, SE╲, SW╱)
# My coordinates: ▽(r,j) has top─, and two sides.
# ▽ vertices: (r, j+1), (r+1, j), (r+1, j+1)
#   left side: (r, j+1)-(r+1, j) → goes UP-LEFT → ╲ direction
#   right side: (r, j+1)-(r+1, j+1) → goes UP-RIGHT → ╱ direction
#   top: (r+1, j)-(r+1, j+1) → horizontal ─
#
# SageMath Nabla: (north─, SE╲, SW╱)
# In my coords: (top─, left╲, right╱)
# So the triple is: (top, left╲, right╱)
NABLA_TRIPLES = {
    (0, 0, 0),
    (0, 10, 1),
    (10, 1, 0),
    (1, 0, 10),
    (1, 1, 1),
}

LABELS = [0, 1, 10]  # possible internal edge labels


def solve(lam, mu, nu):
    """Find all KTW puzzles for c^nu_{lam, mu}.
    Convention (Molev): bottom=nu, left=lam, right=mu."""

    n = max(len(lam), len(mu), len(nu), 1)
    m = max((list(lam)+[0])[0], (list(mu)+[0])[0], (list(nu)+[0])[0], 1)
    N = n + m

    bottom = partition_to_binary(nu, n, m)
    left   = partition_to_binary(lam, n, m)
    right  = partition_to_binary(mu, n, m)

    print(f"n={n}, m={m}, N={N}")
    print(f"bottom (ν={nu}): {bottom}")
    print(f"left   (λ={lam}): {left}")
    print(f"right  (μ={mu}): {right}")

    # Fix boundary edge labels
    boundary = {}
    for j in range(N):
        boundary[edge((0,j),(0,j+1))] = bottom[j]
    for r in range(N):
        boundary[edge((r,0),(r+1,0))] = left[r]
    for k in range(N):
        r = N - 1 - k
        boundary[edge((r, N-r), (r+1, N-r-1))] = right[k]

    # Row-by-row solver. Process row r, assign all new edges, filter by triangles.
    states = [dict(boundary)]

    for r in range(N):
        width = N - r
        if width == 0:
            break

        # Identify new (internal) edges at this row
        new_edges = []
        seen = set()
        for j in range(width):
            for e in [edge((r,j),(r+1,j)), edge((r,j+1),(r+1,j))]:
                if e not in boundary and e not in seen:
                    new_edges.append(e)
                    seen.add(e)
        for j in range(width - 1):
            e = edge((r+1,j),(r+1,j+1))
            if e not in boundary and e not in seen:
                new_edges.append(e)
                seen.add(e)

        n_new = len(new_edges)
        print(f"Row {r}: {width}△ + {width-1}▽, {n_new} new edges, {len(states)} states", end="")

        # Enumerate all label assignments for new edges (each from {0, 1, 10})
        # Use itertools.product for 3^n_new possibilities
        from itertools import product as iprod
        new_states = []

        for assignment_vals in iprod(LABELS, repeat=n_new):
            assignment = dict(zip(new_edges, assignment_vals))

            for state in states:
                merged = dict(state)
                merged.update(assignment)

                ok = True
                # Check △ triangles: triple = (bottom─, left╱, right╲)
                for j in range(width):
                    b = merged[edge((r,j),(r,j+1))]
                    l = merged[edge((r,j),(r+1,j))]
                    rv = merged[edge((r,j+1),(r+1,j))]
                    if (b, l, rv) not in DELTA_TRIPLES:
                        ok = False
                        break
                if not ok:
                    continue

                # Check ▽ triangles: triple = (top─, left╲, right╱)
                for j in range(width - 1):
                    t = merged[edge((r+1,j),(r+1,j+1))]
                    lv = merged[edge((r,j+1),(r+1,j))]
                    rvv = merged[edge((r,j+1),(r+1,j+1))]
                    if (t, rvv, lv) not in NABLA_TRIPLES:
                        ok = False
                        break
                if ok:
                    new_states.append(dict(merged))

        states = new_states
        print(f" → {len(states)}")

        if not states:
            print("  No valid states!")
            break

    print(f"\nTotal solutions: {len(states)}")
    return states, N


def edge_direction(e):
    """Classify edge as 'H' (horizontal), 'slash' (╱), or 'backslash' (╲)."""
    (r1, j1), (r2, j2) = e
    if r1 == r2:
        return 'H'
    # Sorted: lower row first. ╱: same j. ╲: lower row has higher j.
    if j1 == j2:
        return 'slash'     # ╱: (r,j)→(r+1,j)
    else:
        return 'backslash'  # ╲: (r,j+1)→(r+1,j)


def adjacent_edges_at(vertex, N):
    """Return all edges in the triangular grid adjacent to vertex (R, J)."""
    R, J = vertex
    result = []
    # ─ right: (R,J)-(R,J+1)
    if J < N - R:
        result.append(edge((R, J), (R, J + 1)))
    # ─ left: (R,J-1)-(R,J)
    if J > 0:
        result.append(edge((R, J - 1), (R, J)))
    # ╱ up: (R,J)-(R+1,J)
    if R < N and J <= N - R - 1:
        result.append(edge((R, J), (R + 1, J)))
    # ╱ down: (R-1,J)-(R,J)
    if R > 0:
        result.append(edge((R - 1, J), (R, J)))
    # ╲ up-left: (R,J)-(R+1,J-1)
    if R < N and J >= 1 and J <= N - R:
        result.append(edge((R, J), (R + 1, J - 1)))
    # ╲ down-right: (R-1,J+1)-(R,J)
    if R >= 1 and J + 1 <= N - (R - 1):
        result.append(edge((R - 1, J + 1), (R, J)))
    return result


def trace_path(sol, N, start_edge, target_label):
    """Trace a path by following edges with target_label from start_edge.

    At each vertex, looks for exactly one continuation edge with the same
    label (excluding the edge we arrived on). Returns the longest trace
    found when trying both directions from the start edge.
    """
    v1, v2 = start_edge
    best_path = [start_edge]

    for start_v in [v1, v2]:
        curr_v = start_v
        prev_e = start_edge
        trial_path = [start_edge]
        visited = {start_edge}

        while True:
            adj = adjacent_edges_at(curr_v, N)
            nexts = [e for e in adj
                     if e not in visited and e in sol and sol[e] == target_label]
            if len(nexts) == 1:
                ne = nexts[0]
                trial_path.append(ne)
                visited.add(ne)
                nv1, nv2 = ne
                curr_v = nv2 if nv1 == curr_v else nv1
            else:
                # 0 continuations = boundary reached; >1 = ambiguity (shouldn't happen)
                break

        if len(trial_path) > len(best_path):
            best_path = trial_path

    return best_path


def path_to_vertex_coords(path_edges, N):
    """Convert path edges to a list of vertex (x,y) coordinates."""
    h = math.sqrt(3) / 2
    def pos(r, j):
        return (j + r * 0.5, r * h)

    if not path_edges:
        return []
    if len(path_edges) == 1:
        v1, v2 = path_edges[0]
        return [pos(*v1), pos(*v2)]

    # Find start vertex (not shared with second edge)
    e0, e1 = path_edges[0], path_edges[1]
    shared = set(e0) & set(e1)
    start_v = [v for v in e0 if v not in shared]
    if not start_v:
        return []
    start_v = start_v[0]

    coords = [pos(*start_v)]
    curr_v = start_v
    for e in path_edges:
        v1, v2 = e
        next_v = v2 if v1 == curr_v else v1
        coords.append(pos(*next_v))
        curr_v = next_v
    return coords


def generate_tikz(sol, N, title="", mode="labels"):
    """Generate TikZ code for a puzzle solution.

    mode='labels': show edge labels (0, 1, 10) with colors
    mode='paths': show green/red path model
    mode='both': grid + paths overlaid
    """
    h = math.sqrt(3) / 2

    def pos(r, j):
        return (j + r * 0.5, r * h)

    lines = []
    lines.append(r"\begin{tikzpicture}[scale=2.0]")
    if title:
        lines.append(f"  \\node[above] at ({N*0.5:.2f}, {N*h + 0.3:.2f}) {{\\small {title}}};")

    if mode in ('labels', 'both'):
        # Draw all edges with label coloring
        drawn = set()
        for r in range(N):
            for j in range(N - r):
                for v1, v2 in [((r,j),(r,j+1)), ((r,j),(r+1,j)), ((r,j+1),(r+1,j))]:
                    e = edge(v1, v2)
                    if e in drawn:
                        continue
                    drawn.add(e)
                    label = sol[e]
                    x1, y1 = pos(*v1)
                    x2, y2 = pos(*v2)
                    if mode == 'both':
                        # Faint grid
                        lines.append(f"  \\draw[black!15, thin] ({x1:.3f},{y1:.3f}) -- ({x2:.3f},{y2:.3f});")
                    else:
                        if label == 1:
                            color, lw = "blue!80!black", "very thick"
                        elif label == 10:
                            color, lw = "red!70!black", "thick, dashed"
                        else:
                            color, lw = "black!30", "thin"
                        lines.append(f"  \\draw[{color}, {lw}] ({x1:.3f},{y1:.3f}) -- ({x2:.3f},{y2:.3f});")

        if mode == 'labels':
            # Edge labels
            drawn = set()
            for r in range(N):
                for j in range(N - r):
                    for v1, v2, anchor in [
                        ((r,j),(r,j+1), "below"),
                        ((r,j),(r+1,j), "left"),
                        ((r,j+1),(r+1,j), "right"),
                    ]:
                        e = edge(v1, v2)
                        if e in drawn:
                            continue
                        drawn.add(e)
                        label = sol[e]
                        x1, y1 = pos(*v1)
                        x2, y2 = pos(*v2)
                        mx, my = (x1+x2)/2, (y1+y2)/2
                        if label == 1:
                            color = "blue!80!black"
                        elif label == 10:
                            color = "red!70!black"
                        else:
                            color = "black!50"
                        lbl_str = "10" if label == 10 else str(label)
                        lines.append(f"  \\node[{anchor}, font=\\tiny, {color}] at ({mx:.3f},{my:.3f}) {{${lbl_str}$}};")

    if mode in ('paths', 'both'):
        # Draw faint grid as background
        drawn = set()
        for r in range(N):
            for j in range(N - r):
                for v1, v2 in [((r,j),(r,j+1)), ((r,j),(r+1,j)), ((r,j+1),(r+1,j))]:
                    e = edge(v1, v2)
                    if e in drawn:
                        continue
                    drawn.add(e)
                    x1, y1 = pos(*v1)
                    x2, y2 = pos(*v2)
                    lines.append(f"  \\draw[black!12, thin] ({x1:.3f},{y1:.3f}) -- ({x2:.3f},{y2:.3f});")

        # Collect all boundary edges and trace paths from them.
        # Green (label 0): trace from every boundary edge labeled 0.
        # Red   (label 1): trace from every boundary edge labeled 1.
        traced_green = set()
        traced_red = set()
        green_paths = []
        red_paths = []

        # All boundary edges: left ╱, right ╲, bottom ─
        boundary_edges = []
        for r in range(N):
            boundary_edges.append(edge((r, 0), (r + 1, 0)))              # left ╱
        for k in range(N):
            r = N - 1 - k
            boundary_edges.append(edge((r, N - r), (r + 1, N - r - 1)))  # right ╲
        for j in range(N):
            boundary_edges.append(edge((0, j), (0, j + 1)))              # bottom ─

        for be in boundary_edges:
            lbl = sol.get(be, -1)
            if lbl == 0 and be not in traced_green:
                gp = trace_path(sol, N, be, target_label=0)
                for e in gp:
                    traced_green.add(e)
                if len(gp) >= 1:
                    green_paths.append(gp)
            elif lbl == 1 and be not in traced_red:
                rp = trace_path(sol, N, be, target_label=1)
                for e in rp:
                    traced_red.add(e)
                if len(rp) >= 1:
                    red_paths.append(rp)

        off = 0.04  # perpendicular offset so green/red don't overlap

        for path in green_paths:
            coords = path_to_vertex_coords(path, N)
            if len(coords) < 2:
                continue
            segs = " -- ".join(f"({x - off:.3f},{y + off:.3f})" for x, y in coords)
            lines.append(f"  \\draw[green!60!black, line width=2pt, line cap=round, line join=round] {segs};")

        for path in red_paths:
            coords = path_to_vertex_coords(path, N)
            if len(coords) < 2:
                continue
            segs = " -- ".join(f"({x + off:.3f},{y - off:.3f})" for x, y in coords)
            lines.append(f"  \\draw[red!70!black, line width=2pt, line cap=round, line join=round] {segs};")
    # Vertices
    for r in range(N + 1):
        for j in range(N - r + 1):
            x, y = pos(r, j)
            lines.append(f"  \\fill ({x:.3f},{y:.3f}) circle (0.6pt);")

    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines)


if __name__ == "__main__":
    print("=" * 60)
    print("KTW Puzzle Solver (3-label: 0, 1, 10)")
    print("=" * 60)

    examples = [
        ("Molev", (2,), (1,), (2,1), 1),
        ("ZJ", (1,1), (1,1), (2,1,1), 1),
        ("c=2", (2,1), (2,1), (3,2,1), 2),
    ]

    all_results = {}
    for name, lam, mu, nu, expected in examples:
        print(f"\n--- {name}: c^{nu}_{{{lam},{mu}}} = {expected} ---")
        sols, N = solve(lam=lam, mu=mu, nu=nu)
        status = "✓" if len(sols) == expected else "✗"
        print(f"Expected: {expected}, Got: {len(sols)} {status}")
        all_results[name] = (sols, N, lam, mu, nu)

    # Generate LaTeX document: labels row + path model row
    doc = [
        r"\documentclass[border=10pt]{standalone}",
        r"\usepackage{tikz}",
        r"\usepackage{xcolor}",
        r"\usepackage{amsmath}",
        r"\begin{document}",
    ]

    def make_title(name, lam, mu, nu, sols, i=None):
        lam_s = ",".join(map(str, lam))
        mu_s = ",".join(map(str, mu))
        nu_s = ",".join(map(str, nu))
        base = f"$c^{{({nu_s})}}__{{({lam_s}),({mu_s})}}={len(sols)}$"
        if i is not None and len(sols) > 1:
            return f"Puzzle {i+1}/{len(sols)}: {base}"
        return base

    # Row 1: edge-label view
    for name, (sols, N, lam, mu, nu) in all_results.items():
        for i, sol in enumerate(sols):
            title = make_title(name, lam, mu, nu, sols, i if len(sols) > 1 else None)
            doc.append(generate_tikz(sol, N, title=title, mode="labels"))
            doc.append(r"\hspace{0.8cm}")

    doc.append(r"\vspace{1.5cm}")
    doc.append("")

    # Row 2: path model view
    for name, (sols, N, lam, mu, nu) in all_results.items():
        for i, sol in enumerate(sols):
            title = make_title(name, lam, mu, nu, sols, i if len(sols) > 1 else None)
            title = title.replace("$c^", "$\\text{paths: }c^")
            doc.append(generate_tikz(sol, N, title=title, mode="paths"))
            doc.append(r"\hspace{0.8cm}")

    doc.append(r"\end{document}")

    path = "/Users/robin/git/puzzles/puzzles.tex"
    with open(path, 'w') as f:
        f.write("\n".join(doc))
    print(f"\nLaTeX saved to {path}")
