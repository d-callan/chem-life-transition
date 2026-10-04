"""Assembly-index observables C(lambda).

Two variants:

- path_assembly: cheap proxy — shortest joining-path depth from food
  species to each molecule in THIS network (1 + max over reactants).
  Well-defined for any ReactionNetwork. Catalysts are ignored (they are
  not consumed), so an autocatalytic loop does not inflate its own index.
- true assembly index (assembly theory): minimum joining steps over
  molecular substructure. Only defined for generators that attach real
  structure (e.g. polymer_ligation strings). TODO milestone 4.

The key observable is the DISTRIBUTION P(A), not the mean: the
conversation predicts the high-assembly tail gets populated and becomes
persistent only above the transition.
"""

import numpy as np

from ..networks.core import ReactionNetwork


def path_assembly_index(net: ReactionNetwork) -> np.ndarray:
    """Shortest-path assembly index for every species (inf if unreachable).

    A(food) = 0; A(p) = min over producing reactions of
    1 + max(A(reactants)). Relaxation to fixed point — converges because
    estimates only decrease and are bounded below.
    """
    A = np.full(net.n_species, np.inf)
    A[list(net.food)] = 0.0
    for _ in range(10 * net.n_species):
        changed = False
        for rxn in net.reactions:
            src = [A[i] for i in rxn.reactants]
            if not np.all(np.isfinite(src)):
                continue
            cand = 1.0 + max(src)
            for p in rxn.products:
                if cand < A[p]:
                    A[p] = cand
                    changed = True
        if not changed:
            break
    return A


def _string_assembly(s: str, memo: dict) -> int:
    """Assembly index of a string (Marshall/Cronin-style).

    a(s) = min over bipartitions s=u+v of
        a(u) + 1            if u == v   (build once, join to itself)
        a(u) + a(v) + 1     otherwise
    with a(single letter) = 0. Previously built strings stay available,
    so repeats are cheap — that's the whole point of the measure.
    """
    if len(s) <= 1:
        return 0
    if s in memo:
        return memo[s]
    best = len(s) - 1  # monomer-by-monomer upper bound
    for cut in range(1, len(s)):
        u, v = s[:cut], s[cut:]
        cost = _string_assembly(u, memo) + 1 if u == v else \
            _string_assembly(u, memo) + _string_assembly(v, memo) + 1
        best = min(best, cost)
    memo[s] = best
    return best


def molecular_assembly_index(net: ReactionNetwork) -> np.ndarray:
    """True assembly index from species structure (requires structure field)."""
    if not net.structure:
        raise ValueError("no species structure — use path_assembly_index "
                         "or a structured generator (polymer_ligation)")
    memo: dict = {}
    A = np.full(net.n_species, np.inf)
    for i in range(net.n_species):
        s = net.structure.get(i)
        if s is not None:
            A[i] = _string_assembly(s, memo)
    return A


def assembly_distribution(
    net: ReactionNetwork, abundance: np.ndarray, true_index: bool = False
) -> tuple[np.ndarray, np.ndarray]:
    """Abundance-weighted P(A) over producible species.

    Returns (bin_edges, probabilities). Species unreachable from food
    (A=inf) are excluded — they can't be made anyway.
    """
    A = (molecular_assembly_index(net) if true_index
         else path_assembly_index(net))
    abundance = np.asarray(abundance, dtype=float)
    mask = np.isfinite(A)
    if not np.any(mask):
        return np.array([0.0, 1.0]), np.array([1.0])
    A_f, w = A[mask], abundance[mask]
    if w.sum() <= 0:
        w = np.ones_like(w)
    lo, hi = int(A_f.min()), int(A_f.max())
    bins = np.arange(lo, hi + 2) - 0.5
    hist, edges = np.histogram(A_f, bins=bins, weights=w, density=True)
    return edges, hist
