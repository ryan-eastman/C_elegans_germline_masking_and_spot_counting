"""One-to-one minimum-cost matching on a sparse candidate graph.

Both benchmarks match two sets of objects one to one: our spots to the expert's spots within a small
cylinder, our nuclei to the traced nuclei whose masks overlap them. The exact answer is the Hungarian
algorithm (scipy.optimize.linear_sum_assignment), but a dense cost matrix over 30000 by 30000 spots
would take 7 GB. The candidate graph is sparse, so it is split into its connected components first:
no candidate pair joins two components, so solving every component on its own small dense matrix is
exactly the global optimum (a component is usually one or two spots on each side). ``max_dense_cells``
bounds the largest dense matrix; a component above it is matched greedily (nearest pair first) and the
result records that it happened, so a reader knows the score for that gonad is not the exact optimum.
A caller can also cap the total number of candidate pairs (``max_candidates``) and fall back to greedy
matching for the whole set.

Costs must be finite and non-negative. Where several candidate pairs share a cost the result is still
deterministic: ties are broken by the sorted order of the indices, and the dense matrices are built
from sorted unique ids.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

DEFAULT_MAX_DENSE_CELLS = 4_000_000     # 32 MB of float64 per component matrix
DEFAULT_MAX_CANDIDATES = 5_000_000      # candidate pairs per gonad before the greedy fallback


@dataclass
class Assignment:
    """Matched pairs: ``rows[k]`` (first set) is matched with ``cols[k]`` (second set) at ``cost[k]``."""

    rows: np.ndarray
    cols: np.ndarray
    cost: np.ndarray
    n_candidates: int = 0
    n_components: int = 0
    n_greedy_components: int = 0
    method: str = "none"
    notes: list = field(default_factory=list)

    @property
    def n_matched(self) -> int:
        return int(len(self.rows))

    @property
    def total_cost(self) -> float:
        return float(np.sum(self.cost)) if len(self.cost) else 0.0

    def as_frame(self) -> pd.DataFrame:
        return pd.DataFrame({"row": self.rows, "col": self.cols, "cost": self.cost})


def _empty(method: str = "none", n_candidates: int = 0, notes: list | None = None) -> Assignment:
    z = np.zeros(0, dtype=np.int64)
    return Assignment(z, z.copy(), np.zeros(0, float), n_candidates, 0, 0, method, notes or [])


def _dedupe(rows, cols, cost):
    """One entry per (row, col) pair, keeping the smallest cost."""
    order = np.lexsort((cost, cols, rows))
    rows, cols, cost = rows[order], cols[order], cost[order]
    keep = np.ones(len(rows), bool)
    keep[1:] = (rows[1:] != rows[:-1]) | (cols[1:] != cols[:-1])
    return rows[keep], cols[keep], cost[keep]


def big_cost(costs, n_rows: int, n_cols: int) -> float:
    """A finite cost for a forbidden cell that no assignment prefers over one more real match: larger
    than the sum of every real cost an assignment can contain."""
    cmax = float(np.max(costs)) if np.size(costs) else 0.0
    return (min(n_rows, n_cols) + 1.0) * (cmax + 1.0)


def assign_greedy(rows, cols, cost, n_rows: int, n_cols: int) -> np.ndarray:
    """Greedy nearest matching: take candidate pairs in order of cost (ties by row, then column) while
    both ends are free. Returns the indices of the kept candidate pairs."""
    rows = np.asarray(rows, np.int64)
    cols = np.asarray(cols, np.int64)
    cost = np.asarray(cost, float)
    order = np.lexsort((cols, rows, cost))
    used_r = np.zeros(int(n_rows), bool)
    used_c = np.zeros(int(n_cols), bool)
    keep = []
    for k in order:
        r, c = rows[k], cols[k]
        if used_r[r] or used_c[c]:
            continue
        used_r[r] = True
        used_c[c] = True
        keep.append(k)
    return np.asarray(keep, dtype=np.int64)


def assign_dense(cost_matrix: np.ndarray, feasible: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Reference implementation: the Hungarian algorithm on a full matrix in which ``feasible`` marks
    the candidate cells. Used by the tests to check the chunked version and fine for small problems."""
    cost_matrix = np.asarray(cost_matrix, float)
    feasible = np.asarray(feasible, bool)
    if cost_matrix.size == 0 or not feasible.any():
        z = np.zeros(0, np.int64)
        return z, z.copy()
    big = big_cost(cost_matrix[feasible], *cost_matrix.shape)
    cm = np.where(feasible, cost_matrix, big)
    r, c = linear_sum_assignment(cm)
    ok = feasible[r, c]
    return r[ok], c[ok]


def assign_sparse(rows, cols, cost, n_rows: int, n_cols: int, *,
                  max_dense_cells: int = DEFAULT_MAX_DENSE_CELLS,
                  max_candidates: int | None = DEFAULT_MAX_CANDIDATES) -> Assignment:
    """Minimum-cost one-to-one matching over the candidate pairs (rows[k], cols[k]) with cost[k].

    The pairs are grouped into connected components and each component is solved exactly with the
    Hungarian algorithm on its own dense matrix (components with a single row or column are solved by
    picking the cheapest pair). A component whose matrix would exceed ``max_dense_cells`` is matched
    greedily instead and counted in ``n_greedy_components``; more candidate pairs than
    ``max_candidates`` switches the whole problem to greedy matching. ``method`` says which happened.
    """
    rows = np.asarray(rows, np.int64).ravel()
    cols = np.asarray(cols, np.int64).ravel()
    cost = np.asarray(cost, float).ravel()
    if not (len(rows) == len(cols) == len(cost)):
        raise ValueError("rows, cols and cost must have the same length")
    if rows.size == 0:
        return _empty()
    if not np.all(np.isfinite(cost)) or np.any(cost < 0):
        raise ValueError("costs must be finite and non-negative")
    if rows.min() < 0 or cols.min() < 0 or rows.max() >= n_rows or cols.max() >= n_cols:
        raise ValueError("candidate indices out of range")

    rows, cols, cost = _dedupe(rows, cols, cost)
    n_cand = int(len(rows))
    if max_candidates is not None and n_cand > max_candidates:
        keep = assign_greedy(rows, cols, cost, n_rows, n_cols)
        note = (f"{n_cand} candidate pairs exceed the budget of {max_candidates}; "
                "greedy nearest matching was used for the whole set")
        return Assignment(rows[keep], cols[keep], cost[keep], n_cand, 0, 0, "greedy", [note])

    n_nodes = int(n_rows) + int(n_cols)
    graph = coo_matrix((np.ones(n_cand), (rows, cols + n_rows)), shape=(n_nodes, n_nodes))
    _, comp = connected_components(graph, directed=False)
    edge_comp = comp[rows]
    order = np.argsort(edge_comp, kind="stable")
    bounds = np.flatnonzero(np.diff(edge_comp[order])) + 1
    groups = np.split(order, bounds)

    selected: list[int] = []
    n_greedy = 0
    for g in groups:
        r, c, w = rows[g], cols[g], cost[g]
        ur, uc = np.unique(r), np.unique(c)
        if len(ur) == 1 or len(uc) == 1:
            selected.append(int(g[np.lexsort((c, r, w))[0]]))
            continue
        if len(ur) * len(uc) > max_dense_cells:
            n_greedy += 1
            keep = assign_greedy(r, c, w, n_rows, n_cols)
            selected.extend(int(k) for k in g[keep])
            continue
        ri = np.searchsorted(ur, r)
        ci = np.searchsorted(uc, c)
        big = big_cost(w, len(ur), len(uc))
        m = np.full((len(ur), len(uc)), big)
        m[ri, ci] = w
        edge_index = np.full(m.shape, -1, np.int64)
        edge_index[ri, ci] = g
        ar, ac = linear_sum_assignment(m)
        ok = edge_index[ar, ac] >= 0
        selected.extend(int(k) for k in edge_index[ar[ok], ac[ok]])

    sel = np.array(sorted(selected), dtype=np.int64)
    method = "hungarian" if n_greedy == 0 else "hungarian+greedy"
    notes = [] if n_greedy == 0 else [
        f"{n_greedy} component(s) exceeded max_dense_cells={max_dense_cells} and were matched greedily"]
    return Assignment(rows[sel], cols[sel], cost[sel], n_cand, len(groups), n_greedy, method, notes)
