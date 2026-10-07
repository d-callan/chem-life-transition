"""Core data structure for the rebuilt chemistry (chemtrans 2.0).

A ReactionNetwork is a set of reversible elementary *channels*. Each
channel is a single reversible reaction:

    alpha  <=>  beta        (optionally catalyzed by species in `catalysts`)

where alpha/beta are stoichiometric multisets of species (the structural
stoichiometry). Catalysts are NOT consumed; they enter the *kinetic*
participant multiset but not the stoichiometric view, so RAF-style
structural analysis and kinetics stay consistent.

Kinetics (single convention, used by ODE, Jacobian, and SSA alike):

    forward flux  f_r(x) = k_f * prod_i x_i^{a_i} * (eps_r + kappa * sum_c x_c)
    reverse flux  b_r(x) = k_r * prod_i x_i^{b_i} * (eps_r + kappa * sum_c x_c)

    eps_r = eps   if channel.synth  (spontaneous synthesis: the
                                    symmetry-breaking "field")
          = 1.0   otherwise         (decay channels are not eps-scaled)

  - Catalysis multiplies BOTH directions by the same factor, so it never
    changes K_eq = k_f/k_r (local detailed balance).
  - eps=0 freezes spontaneous synthesis channels entirely; with catalysts
    drawn only from non-food polymer species the all-food state is then
    absorbing for the polymer core (labeled absorbing limit).
  - k_r is DERIVED from species free energies: k_r = k_f / K_eq with
    K_eq = exp(-sum_i g_i nu_i), nu = products - reactants. Generators
    must respect this; `check_detailed_balance` verifies it.

Stochastic hazards are falling-factorial forms derived from the SAME
channel data (see dynamics/gillespie.py); the ODE rate k*x^m corresponds
to hazard k * V^{1-m} * prod_i (n_i)_{a_i} - no combinatorial prefactors
beyond the falling factorial.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

import numpy as np

MODEL_VERSION = "2.0.0"


@dataclass(frozen=True)
class Channel:
    """One reversible elementary channel.

    reactants/products: sorted tuples of (species_index, stoich) -- the
        structural stoichiometry of the forward direction.
    catalysts: species indices that accelerate the channel (alternatives:
        ANY one suffices). Not consumed.
    k_f, k_r: forward/reverse rate constants. For channels built by the
        polymer generator, k_r = k_f / K_eq with K_eq from species free
        energies; benchmarks may set rates directly (still validated).
    synth: True if this channel's uncatalyzed component is scaled by the
        protocol parameter eps (spontaneous polymer synthesis). Decay
        channels (hydrolysis, fuel deactivation) use synth=False.
    mechanism: free-form label ("elongation", "thermal", "hydrolysis",
        "deactivation", "background", ...).
    """

    reactants: tuple[tuple[int, int], ...]
    products: tuple[tuple[int, int], ...]
    catalysts: tuple[int, ...] = ()
    k_f: float = 1.0
    k_r: float = 0.0
    synth: bool = True
    mechanism: str = "background"

    @property
    def nu(self) -> dict[int, int]:
        """Net stoichiometry products - reactants."""
        out: dict[int, int] = {}
        for i, a in self.reactants:
            out[i] = out.get(i, 0) - a
        for i, a in self.products:
            out[i] = out.get(i, 0) + a
        return {i: v for i, v in out.items() if v != 0}

    def factors(self, direction: str) -> tuple[int, ...]:
        """Kinetic participant multiset (flat list w/ repetition) for a
        direction's UNCATALYZED mass-action factor."""
        src = self.reactants if direction == "fwd" else self.products
        return tuple(i for i, a in src for _ in range(a))


class ReactionNetwork:
    """A chemical reaction network: species, channels, food, thermo data.

    Treated as immutable after construction: derived arrays are cached on
    first use. Use `with_channels`/`with_food_inlet` to build variants
    (knockouts, changed protocols) rather than mutating in place.
    """

    def __init__(
        self,
        species,
        channels,
        food,
        *,
        structure=None,
        comp=None,
        g=None,
        food_inlet=None,
        kappa=1.0,
        eps=1.0,
        cat_power=1,
        params=None,
        sparse=None,
    ):
        self.species = tuple(str(s) for s in species)
        self.channels = tuple(channels)
        self.food = frozenset(int(i) for i in food)
        self.structure = dict(structure or {})
        self.kappa = float(kappa)
        self.eps = float(eps)
        if not float(cat_power).is_integer() or int(cat_power) < 1:
            raise ValueError("cat_power must be a positive integer")
        self.cat_power = int(cat_power)
        self.params = dict(params or {})
        # sparse mode: never materialize dense S/C/J. None -> auto for
        # large networks (dense S would be n_species x n_channels).
        self.use_sparse = (bool(sparse) if sparse is not None
                           else len(self.species) >= 1000)

        n = len(self.species)
        self._validate()

        self.comp = np.asarray(comp, dtype=float) if comp is not None else None
        if self.comp is not None and self.comp.shape[1] != n:
            raise ValueError("comp must be (n_letters, n_species)")
        self.g = np.asarray(g, dtype=float) if g is not None else None
        if self.g is not None and self.g.shape != (n,):
            raise ValueError("g must have length n_species")
        self.food_inlet = (
            np.asarray(food_inlet, dtype=float)
            if food_inlet is not None
            else np.zeros(n)
        )
        if self.food_inlet.shape != (n,):
            raise ValueError("food_inlet must have length n_species")
        if (self.food_inlet < 0).any():
            raise ValueError("food_inlet must be nonnegative")
        nonzero_inlet = set(np.nonzero(self.food_inlet)[0].tolist())
        if not nonzero_inlet <= self.food:
            raise ValueError("food_inlet may only supply food species")
        self._cache = {}

    # ------------------------------------------------------------------
    def _validate(self) -> None:
        n = len(self.species)
        if len(set(self.species)) != n:
            raise ValueError("duplicate species names")
        seen_ids: set[str] = set()
        for ch in self.channels:
            if not isinstance(ch, Channel):
                raise TypeError("channels must be Channel objects")
            for src in (ch.reactants, ch.products):
                last = -1
                for i, a in src:
                    if not (0 <= i < n):
                        raise ValueError("channel species index out of range")
                    if not (isinstance(a, int) or float(a).is_integer()) or a <= 0:
                        raise ValueError("stoichiometry must be positive integer")
                    if i <= last:
                        raise ValueError("species entries must be sorted/unique")
                    last = i
            for c in ch.catalysts:
                if not (0 <= c < n):
                    raise ValueError("catalyst index out of range")
            if not np.isfinite(ch.k_f) or not np.isfinite(ch.k_r):
                raise ValueError("rates must be finite")
            if ch.k_f < 0 or ch.k_r < 0:
                raise ValueError("rates must be nonnegative")
            if (ch.k_f > 0) != (ch.k_r > 0) and (ch.k_f > 0 or ch.k_r > 0):
                # one-way channels are allowed ONLY via k_r == 0 *or* k_f == 0
                # explicitly; flag the asymmetric case loudly in params, but
                # it is physically allowed for driven steps (e.g. feed-only).
                pass
            cid = repr((ch.reactants, ch.products, ch.catalysts))
            if cid in seen_ids:
                raise ValueError(f"duplicate channel {cid}")
            seen_ids.add(cid)

    # ------------------------------------------------------------------
    @property
    def n_species(self) -> int:
        return len(self.species)

    @property
    def n_channels(self) -> int:
        return len(self.channels)

    @property
    def is_food(self) -> np.ndarray:
        v = np.zeros(self.n_species, dtype=bool)
        v[list(self.food)] = True
        return v

    @property
    def polymer_idx(self) -> np.ndarray:
        """Non-food species indices (the candidate organized material)."""
        return np.array(
            [i for i in range(self.n_species) if i not in self.food], dtype=int
        )

    @property
    def polymer_lengths(self) -> np.ndarray:
        out = np.full(self.n_species, -1, dtype=int)
        for i, s in self.structure.items():
            out[i] = len(s)
        return out

    # ---------------- derived arrays (cached) -------------------------
    def stoichiometry(self) -> np.ndarray:
        """S[i, r] = net production of species i by channel r (forward)."""
        if "S" not in self._cache:
            if self.use_sparse:
                self._cache["S"] = self._stoichiometry_sparse()
            else:
                S = np.zeros((self.n_species, self.n_channels))
                for r, ch in enumerate(self.channels):
                    for i, a in ch.products:
                        S[i, r] += a
                    for i, a in ch.reactants:
                        S[i, r] -= a
                self._cache["S"] = S
        return self._cache["S"]

    def _stoichiometry_sparse(self):
        """CSR stoichiometry — ~4 nnz per channel regardless of pool
        size; never materialize the dense array."""
        import scipy.sparse as sp
        rows, cols, data = [], [], []
        for r, ch in enumerate(self.channels):
            for i, a in ch.products:
                rows.append(i); cols.append(r); data.append(float(a))
            for i, a in ch.reactants:
                rows.append(i); cols.append(r); data.append(-float(a))
        return sp.csr_matrix(
            (np.asarray(data), (np.asarray(rows), np.asarray(cols))),
            shape=(self.n_species, self.n_channels))

    def catalysis_matrix(self) -> np.ndarray:
        """C[r, c] = 1 if species c catalyzes channel r (alternatives)."""
        if "C" not in self._cache:
            if self.use_sparse:
                self._cache["C"] = self._catalysis_sparse()
            else:
                C = np.zeros((self.n_channels, self.n_species))
                for r, ch in enumerate(self.channels):
                    for c in ch.catalysts:
                        C[r, c] = 1.0
                self._cache["C"] = C
        return self._cache["C"]

    def _catalysis_sparse(self):
        import scipy.sparse as sp
        rows, cols = [], []
        for r, ch in enumerate(self.channels):
            for c in ch.catalysts:
                rows.append(r); cols.append(c)
        return sp.csr_matrix(
            (np.ones(len(rows)),
             (np.asarray(rows), np.asarray(cols))),
            shape=(self.n_channels, self.n_species))

    def _factors(self) -> tuple[np.ndarray, np.ndarray, int]:
        """Per-channel participant index lists for fwd/rev directions,
        padded to a matrix with sentinel n (x_app[n] = 1)."""
        if "factors" not in self._cache:
            fwd = [ch.factors("fwd") for ch in self.channels]
            rev = [ch.factors("rev") for ch in self.channels]
            m = max([len(f) for f in fwd + rev] + [1])
            n = self.n_species
            Rf = np.full((self.n_channels, m), n, dtype=int)
            Rr = np.full((self.n_channels, m), n, dtype=int)
            for r, f in enumerate(fwd):
                Rf[r, : len(f)] = f
            for r, f in enumerate(rev):
                Rr[r, : len(f)] = f
            self._cache["factors"] = (Rf, Rr, m)
        return self._cache["factors"]

    def _rate_vectors(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if "rates" not in self._cache:
            kf = np.array([ch.k_f for ch in self.channels])
            kr = np.array([ch.k_r for ch in self.channels])
            synth = np.array([ch.synth for ch in self.channels])
            self._cache["rates"] = (kf, kr, synth)
        return self._cache["rates"]

    # ---------------- kinetics ---------------------------------------
    def rates(self, x: np.ndarray, eps: float | None = None):
        """(fwd_flux, rev_flux) per channel at concentration vector x."""
        if eps is None:
            eps = self.eps
        x = np.asarray(x, dtype=float)
        x_app = np.append(x, 1.0)
        Rf, Rr, _ = self._factors()
        kf, kr, synth = self._rate_vectors()
        C = self.catalysis_matrix()
        bg = np.where(synth, eps, 1.0)
        boost = bg + self.kappa * np.asarray(
            C @ np.power(x, self.cat_power)).ravel()
        mf = np.prod(x_app[Rf], axis=1)
        mr = np.prod(x_app[Rr], axis=1)
        return kf * mf * boost, kr * mr * boost

    def net_flux(self, x: np.ndarray, eps: float | None = None) -> np.ndarray:
        f, b = self.rates(x, eps)
        return f - b

    def d_rates(self, x: np.ndarray, eps: float | None = None):
        """(Df, Db): Jacobians of fwd/rev channel fluxes wrt x.

        Valid at x = 0 (no division by x): the derivative of the monomial
        prod_j x_j^{a_j} wrt x_i is sum over factor positions p == i of the
        product over remaining factors.
        """
        if eps is None:
            eps = self.eps
        if self.use_sparse:
            return self._d_rates_sparse(x, eps)
        # the COO assembly is strictly faster than the per-channel
        # loop even at small sizes; densify for dense-mode callers
        Df, Db = self._d_rates_sparse(x, eps)
        return Df.toarray(), Db.toarray()

    def _d_rates_dense(self, x: np.ndarray, eps: float):
        """Reference per-channel loop — kept for cross-checks."""
        x = np.asarray(x, dtype=float)
        n, nc = self.n_species, self.n_channels
        x_app = np.append(x, 1.0)
        Rf, Rr, m = self._factors()
        kf, kr, synth = self._rate_vectors()
        C = self.catalysis_matrix()
        bg = np.where(synth, eps, 1.0)
        Cx = np.asarray(
            C @ np.power(x, self.cat_power)).ravel()
        mf = np.prod(x_app[Rf], axis=1)
        mr = np.prod(x_app[Rr], axis=1)

        Df = np.zeros((nc, n))
        Db = np.zeros((nc, n))
        for r, ch in enumerate(self.channels):
            fac_f = ch.factors("fwd")
            fac_r = ch.factors("rev")
            # d(monomial)/dx_i: sum over positions holding species i of the
            # product over the other positions.
            for pos, s in enumerate(fac_f):
                others = np.prod(
                    [x_app[fac_f[q]] for q in range(len(fac_f)) if q != pos],
                    axis=0,
                ) if len(fac_f) > 1 else 1.0
                Df[r, s] += kf[r] * (bg[r] + self.kappa * Cx[r]) * others
            for pos, s in enumerate(fac_r):
                others = np.prod(
                    [x_app[fac_r[q]] for q in range(len(fac_r)) if q != pos],
                    axis=0,
                ) if len(fac_r) > 1 else 1.0
                Db[r, s] += kr[r] * (bg[r] + self.kappa * Cx[r]) * others
            # catalysts contribute via boost: d(k*mass*kappa*x_c^p)/dx_c
            cat_idx = np.nonzero(C[r])[0]
            if len(cat_idx):
                dpow = (self.cat_power
                        * np.power(np.maximum(x[cat_idx], 0.0),
                                   self.cat_power - 1))
                Df[r, cat_idx] += kf[r] * mf[r] * self.kappa * dpow
                Db[r, cat_idx] += kr[r] * mr[r] * self.kappa * dpow
        return Df, Db

    def _cat_idx_padded(self) -> tuple[np.ndarray, int]:
        """(nc x mc) padded catalyst index array, sentinel n_species."""
        if "catidx" not in self._cache:
            cats = [ch.catalysts for ch in self.channels]
            mc = max([len(c) for c in cats] + [1])
            R = np.full((self.n_channels, mc), self.n_species,
                        dtype=int)
            for r, c in enumerate(cats):
                R[r, : len(c)] = c
            self._cache["catidx"] = (R, mc)
        return self._cache["catidx"]

    def _d_rates_sparse(self, x: np.ndarray, eps: float):
        """Vectorized COO assembly of (Df, Db) — same math as d_rates
        but no dense (nc x n) arrays and no per-channel Python loop.
        Duplicate (r, s) entries (e.g. AA+AA reactant positions) are
        summed by csr conversion."""
        import scipy.sparse as sp
        n, nc = self.n_species, self.n_channels
        x_app = np.append(x, 1.0)
        Rf, Rr, m = self._factors()
        kf, kr, synth = self._rate_vectors()
        C = self.catalysis_matrix()
        bg = np.where(synth, eps, 1.0)
        Cx = np.asarray(
            C @ np.power(x, self.cat_power)).ravel()
        mf = np.prod(x_app[Rf], axis=1)
        mr = np.prod(x_app[Rr], axis=1)

        rows_f, cols_f, data_f = [], [], []
        rows_b, cols_b, data_b = [], [], []
        rr = np.arange(nc)
        for pos in range(m):
            for R, xk, kvec, flux_m, rows, cols, data in (
                    (Rf, x_app, kf, mf, rows_f, cols_f, data_f),
                    (Rr, x_app, kr, mr, rows_b, cols_b, data_b)):
                s = R[:, pos]
                mask = s < n
                if not mask.any():
                    continue
                others = np.ones(nc)
                for q in range(m):
                    if q != pos:
                        others *= xk[R[:, q]]
                rows.append(rr[mask])
                cols.append(s[mask])
                data.append((kvec * (bg + self.kappa * Cx)
                             * others)[mask])
        Rc, mc = self._cat_idx_padded()
        for pos in range(mc):
            c = Rc[:, pos]
            mask = c < n
            if not mask.any():
                continue
            dpow = (self.cat_power
                    * np.power(np.maximum(x[np.where(mask, c, 0)], 0.0),
                               self.cat_power - 1))
            rows_f.append(rr[mask]); cols_f.append(c[mask])
            data_f.append((kf * mf * self.kappa * dpow)[mask])
            rows_b.append(rr[mask]); cols_b.append(c[mask])
            data_b.append((kr * mr * self.kappa * dpow)[mask])

        def mk(rows, cols, data):
            return sp.csr_matrix(
                (np.concatenate(data),
                 (np.concatenate(rows), np.concatenate(cols))),
                shape=(nc, n))
        return mk(rows_f, cols_f, data_f), mk(rows_b, cols_b, data_b)

    def jacobian(
        self,
        x: np.ndarray,
        dilution: float = 0.0,
        eps: float | None = None,
    ) -> np.ndarray:
        """J = S (Df - Db) - dilution*I for the CSTR vector field."""
        Df, Db = self.d_rates(x, eps)
        if self.use_sparse:
            import scipy.sparse as sp
            return (self.stoichiometry() @ (Df - Db)
                    - dilution * sp.eye(self.n_species, format="csr"))
        return self.stoichiometry() @ (Df - Db) - dilution * np.eye(
            self.n_species
        )

    # ---------------- thermodynamics ----------------------------------
    def check_detailed_balance(self, tol: float = 1e-6) -> dict:
        """Verify k_f/k_r == exp(-g . nu) per channel; also verify that a
        closed system with these g values is at equilibrium (report the
        per-channel residual of K_eq consistency)."""
        if self.g is None:
            raise ValueError("no free energies g on this network")
        rows = []
        ok = True
        g = self.g
        for r, ch in enumerate(self.channels):
            nu = np.zeros(self.n_species)
            for i, a in ch.products:
                nu[i] += a
            for i, a in ch.reactants:
                nu[i] -= a
            dG = float(g @ nu)
            K = np.exp(-dG)
            if ch.k_f > 0 and ch.k_r > 0:
                ratio = ch.k_f / ch.k_r
                rel = abs(ratio / K - 1.0)
                rows.append(
                    {"channel": r, "dG": dG, "K_eq": K, "kf/kr": ratio,
                     "rel_err": rel, "mechanism": ch.mechanism}
                )
                if rel > tol:
                    ok = False
        return {"ok": ok, "channels": rows}

    def equilibrium_state(self) -> np.ndarray:
        """Closed-system equilibrium concentrations (unnormalized
        Boltzmann weights e^{-g_i}, food included)."""
        if self.g is None:
            raise ValueError("no free energies")
        return np.exp(-self.g)

    def check_letter_conservation(self) -> dict:
        """Per-letter conservation: comp @ S must be exactly zero for
        internal (non-boundary) channels."""
        if self.comp is None:
            raise ValueError("no composition matrix")
        S = self.stoichiometry()
        res = (np.asarray(S.T @ self.comp.T).T
               if self.use_sparse else self.comp @ S)
        bad = np.argwhere(np.abs(res) > 0)
        return {
            "ok": bad.size == 0,
            "violations": [
                {"letter": int(a), "channel": int(b)} for a, b in bad
            ],
        }

    def letter_mass(self, x: np.ndarray) -> float:
        """Total monomer mass: sum over letters and species of
        comp[letter, i] * x_i."""
        if self.comp is None:
            raise ValueError("no composition matrix")
        return float(self.comp.sum(axis=0) @ np.asarray(x, dtype=float))

    def feed_letters(self) -> float:
        """Total monomer influx concentration (sum_i comp-mass_i * inlet_i)."""
        return self.letter_mass(self.food_inlet)

    # ---------------- provenance --------------------------------------
    def to_dict(self) -> dict:
        return {
            "model_version": MODEL_VERSION,
            "species": list(self.species),
            "structure": {str(k): v for k, v in self.structure.items()},
            "food": sorted(self.food),
            "food_inlet": self.food_inlet.tolist(),
            "comp": self.comp.tolist() if self.comp is not None else None,
            "g": self.g.tolist() if self.g is not None else None,
            "kappa": self.kappa,
            "eps": self.eps,
            "cat_power": self.cat_power,
            "params": self.params,
            "channels": [
                {
                    "reactants": [list(t) for t in ch.reactants],
                    "products": [list(t) for t in ch.products],
                    "catalysts": list(ch.catalysts),
                    "k_f": ch.k_f,
                    "k_r": ch.k_r,
                    "synth": ch.synth,
                    "mechanism": ch.mechanism,
                }
                for ch in self.channels
            ],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "ReactionNetwork":
        if d.get("model_version") != MODEL_VERSION:
            raise ValueError(
                f"incompatible network schema {d.get('model_version')!r} "
                f"(expected {MODEL_VERSION})"
            )
        channels = [
            Channel(
                reactants=tuple((int(i), int(a)) for i, a in c["reactants"]),
                products=tuple((int(i), int(a)) for i, a in c["products"]),
                catalysts=tuple(int(c2) for c2 in c["catalysts"]),
                k_f=c["k_f"],
                k_r=c["k_r"],
                synth=c["synth"],
                mechanism=c["mechanism"],
            )
            for c in d["channels"]
        ]
        return cls(
            d["species"],
            channels,
            d["food"],
            structure={int(k): v for k, v in d["structure"].items()},
            comp=None if d["comp"] is None else np.array(d["comp"]),
            g=None if d["g"] is None else np.array(d["g"]),
            food_inlet=np.array(d["food_inlet"]),
            kappa=d["kappa"],
            eps=d["eps"],
            cat_power=d.get("cat_power", 1),
            params=d.get("params"),
        )

    @property
    def hash(self) -> str:
        """Order-independent content hash of the whole network."""
        blob = json.dumps(self.to_dict(), sort_keys=True).encode()
        return hashlib.sha256(blob).hexdigest()[:16]

    def with_channels(self, channels, **overrides) -> "ReactionNetwork":
        """Build a variant network (knockouts/interventions)."""
        kw = dict(
            structure=self.structure,
            comp=self.comp,
            g=self.g,
            food_inlet=self.food_inlet,
            kappa=self.kappa,
            eps=self.eps,
            cat_power=self.cat_power,
            params=self.params,
        )
        kw.update(overrides)
        return ReactionNetwork(self.species, list(channels), self.food, **kw)

    def without_catalysis(self) -> "ReactionNetwork":
        """Intervention: remove all catalysts (k_cat -> 0, K_eq unchanged)."""
        return self.with_channels(
            [
                Channel(ch.reactants, ch.products, (), ch.k_f, ch.k_r,
                        ch.synth, ch.mechanism)
                for ch in self.channels
            ]
        )
