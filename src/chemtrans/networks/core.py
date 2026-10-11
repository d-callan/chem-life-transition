"""Core data structure for the rebuilt chemistry (chemtrans 2.0).

A ReactionNetwork is a set of reversible elementary *channels*. Each
channel is a single reversible reaction:

    alpha  <=>  beta        (optionally catalyzed by species in `catalysts`)

where alpha/beta are stoichiometric multisets of species (the structural
stoichiometry). Catalysts are NOT consumed; they enter the *kinetic*
participant multiset but not the stoichiometric view, so RAF-style
structural analysis and kinetics stay consistent.

Kinetics (single convention, used by ODE, Jacobian, and SSA alike):

    forward flux  f_r(x) = k_f * prod_i x_i^{a_i} * (bg_f + kappa * sum_c x_c^p)
    reverse flux  b_r(x) = k_r * prod_i x_i^{b_i} * (bg_r + kappa * sum_c x_c^p)

    bg_d = eps    if channel.synth AND direction d is in
                  channel.synth_dir  (spontaneous synthesis: the
                  symmetry-breaking "field")
         = 1.0    otherwise         (decay directions are not eps-scaled)

  - synth_dir selects WHICH direction carries the eps gate: "both"
    (legacy: synthesis channels gated symmetrically), "fwd", or "rev".
    A channel whose synthesis direction is the REVERSE (e.g. thermal
    PA <=> P + A, where fwd is hydrolysis decay and rev is condensation
    synthesis) uses synth_dir="rev": hydrolysis runs unconditionally at
    eps=0 while condensation stays gated — the dead state remains
    absorbing AND decay channels still run (Invariant 3).

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
    synth_dir: which direction the eps gate applies to when synth=True —
        "both" (default; legacy symmetric gating), "fwd", or "rev".
        Use a directional gate when synthesis and decay live in opposite
        directions of one reversible channel (e.g. thermal
        PA <=> P + A with synth_dir="rev": condensation gated,
        hydrolysis always runs).
    driven: True marks a deliberately irreversible driven-limit channel
        (k_f == 0 or k_r == 0). check_detailed_balance reports it but
        does not fail on it; an undeclared irreversible channel fails.
    mechanism: free-form label ("elongation", "thermal", "hydrolysis",
        "deactivation", "background", ...).
    cat_w: optional per-catalyst edge weights, aligned with `catalysts`
        (empty means every edge weight is 1.0). Used by physically
        derived catalysis — e.g. template affinity weights on ligation
        channels — where catalyst potency depends on the (channel,
        catalyst) pair rather than the catalyst alone. The edge weight
        multiplies kappa*cat_weights[c] identically in both directions,
        so detailed balance is preserved.
    """

    reactants: tuple[tuple[int, int], ...]
    products: tuple[tuple[int, int], ...]
    catalysts: tuple[int, ...] = ()
    k_f: float = 1.0
    k_r: float = 0.0
    synth: bool = True
    mechanism: str = "background"
    synth_dir: str = "both"
    driven: bool = False
    cat_w: tuple[float, ...] = ()

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
        cat_len_power=0.0,
        k_inhibit=0.0,
        inlet_species=None,
        params=None,
        sparse=None,
    ):
        self.species = tuple(str(s) for s in species)
        self.channels = tuple(channels)
        self.food = frozenset(int(i) for i in food)
        # extra species the CSTR inlet may supply beyond food — e.g.
        # environmental short oligomers fed from a bulk background
        # chemistry (f_dimer in the generator). They remain
        # non-food for order-parameter purposes (polymer_mask,
        # is_food unchanged); this only relaxes the inlet check.
        self.inlet_species = frozenset(
            int(i) for i in (inlet_species or ()))
        self.structure = dict(structure or {})
        self.kappa = float(kappa)
        self.eps = float(eps)
        if not float(cat_power).is_integer() or int(cat_power) < 1:
            raise ValueError("cat_power must be a positive integer")
        self.cat_power = int(cat_power)
        if not np.isfinite(cat_len_power) or cat_len_power < 0:
            raise ValueError("cat_len_power must be a nonnegative float")
        self.cat_len_power = float(cat_len_power)
        if not np.isfinite(k_inhibit) or k_inhibit < 0:
            raise ValueError("k_inhibit must be a nonnegative float")
        self.k_inhibit = float(k_inhibit)
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
        allowed_inlet = self.food | self.inlet_species
        if not nonzero_inlet <= allowed_inlet:
            raise ValueError(
                "food_inlet may only supply food or declared "
                "inlet_species")
        self._cache = {}
        # per-species catalytic potency weight (L/2)^cat_len_power for
        # structural polymers (len>=2), 1.0 for everything else — food,
        # activated monomers, waste, and unstructured species. At
        # cat_len_power=0 every weight is exactly 1.0 (legacy flat-kappa
        # behavior). Symmetric across both reaction directions, so
        # detailed balance is untouched.
        self.cat_weights = np.power(
            np.clip(self.polymer_lengths, 2, None) / 2.0,
            self.cat_len_power)

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
            if ch.cat_w:
                if len(ch.cat_w) != len(ch.catalysts):
                    raise ValueError(
                        "cat_w must align with catalysts "
                        f"(got {len(ch.cat_w)} weights for "
                        f"{len(ch.catalysts)} catalysts)")
                if any((not np.isfinite(w)) or w <= 0 for w in ch.cat_w):
                    raise ValueError("cat_w entries must be positive "
                                     "finite floats")
            if not np.isfinite(ch.k_f) or not np.isfinite(ch.k_r):
                raise ValueError("rates must be finite")
            if ch.k_f < 0 or ch.k_r < 0:
                raise ValueError("rates must be nonnegative")
            if ch.synth_dir not in ("fwd", "rev", "both"):
                raise ValueError(
                    f"synth_dir must be 'fwd'|'rev'|'both', "
                    f"got {ch.synth_dir!r}")
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

    @property
    def polymer_mask(self) -> np.ndarray:
        """True for structural polymers (len >= 2); False for food,
        activated monomers, waste, and species without structure
        (e.g. benchmark reservoirs) — the set affected by CSTR
        polymer retention."""
        if "polymer_mask" not in self._cache:
            self._cache["polymer_mask"] = self.polymer_lengths >= 2
        return self._cache["polymer_mask"]

    def dilution_outflow(self, dilution: float, retention: float = 0.0,
                         len_scale: float = 0.0) -> np.ndarray:
        """Per-species CSTR outflow rate vector D_i.

        retention (rho): flat polymer-retention fraction.
        len_scale (gamma in [0,1]): size-selective dialysis mixing —
          rho_i(L) = rho * (1 - 2*gamma/L) for structural polymers.
          gamma=0 -> rho_i = rho for all polymers (legacy flat).
          gamma=1 -> rho_i = rho*(1 - 2/L): dimers wash out at full
          D, longer chains retained progressively (L6 -> 2rho/3).
        Food/fuel/waste and unstructured species always outflow at D.
        """
        if retention:
            if not (0.0 <= len_scale <= 1.0):
                raise ValueError("retention_len_scale must be in [0, 1]")
            lens = np.clip(self.polymer_lengths, 2, None)
            rho_i = (retention * (1.0 - len_scale * 2.0 / lens)
                     * self.polymer_mask)
            return dilution * (1.0 - rho_i)
        return np.full(self.n_species, dilution)

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
        """C[r, c] = catalytic edge weight of species c on channel r.

        1.0 for a plain catalyst; Channel.cat_w supplies a per-edge
        weight when potency depends on the (channel, catalyst) pair
        (e.g. template complementarity on ligation channels)."""
        if "C" not in self._cache:
            if self.use_sparse:
                self._cache["C"] = self._catalysis_sparse()
            else:
                C = np.zeros((self.n_channels, self.n_species))
                for r, ch in enumerate(self.channels):
                    for k, c in enumerate(ch.catalysts):
                        C[r, c] = (ch.cat_w[k] if ch.cat_w else 1.0)
                self._cache["C"] = C
        return self._cache["C"]

    def _catalysis_sparse(self):
        import scipy.sparse as sp
        rows, cols, data = [], [], []
        for r, ch in enumerate(self.channels):
            for k, c in enumerate(ch.catalysts):
                rows.append(r); cols.append(c)
                data.append(ch.cat_w[k] if ch.cat_w else 1.0)
        return sp.csr_matrix(
            (np.asarray(data),
             (np.asarray(rows), np.asarray(cols))),
            shape=(self.n_channels, self.n_species))

    def _catalysis_binary(self):
        """Unweighted 0/1 catalyst incidence — used where catalyst
        *concentrations* (not edge-weighted activity) are needed,
        e.g. the k_inhibit denominator sum_{c in C_r} x_c."""
        if "Cbin" not in self._cache:
            if self.use_sparse:
                import scipy.sparse as sp
                rows, cols = [], []
                for r, ch in enumerate(self.channels):
                    for c in ch.catalysts:
                        rows.append(r); cols.append(c)
                self._cache["Cbin"] = sp.csr_matrix(
                    (np.ones(len(rows)),
                     (np.asarray(rows), np.asarray(cols))),
                    shape=(self.n_channels, self.n_species))
            else:
                Cb = np.zeros((self.n_channels, self.n_species))
                for r, ch in enumerate(self.channels):
                    for c in ch.catalysts:
                        Cb[r, c] = 1.0
                self._cache["Cbin"] = Cb
        return self._cache["Cbin"]

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

    def _rate_vectors(self) -> tuple[np.ndarray, np.ndarray,
                                     np.ndarray, np.ndarray]:
        """(kf, kr, gate_f, gate_r): per-direction eps-gate flags.

        gate_d[r] is True iff channel r's direction-d uncatalyzed
        background is eps-scaled. synth_dir="rev" gates only the
        reverse (e.g. thermal condensation) while the forward decay
        direction runs unconditionally."""
        if "rates" not in self._cache:
            kf = np.array([ch.k_f for ch in self.channels])
            kr = np.array([ch.k_r for ch in self.channels])
            gate_f = np.array([ch.synth and ch.synth_dir != "rev"
                               for ch in self.channels])
            gate_r = np.array([ch.synth and ch.synth_dir != "fwd"
                               for ch in self.channels])
            self._cache["rates"] = (kf, kr, gate_f, gate_r)
        return self._cache["rates"]

    # ---------------- kinetics ---------------------------------------
    def rates(self, x: np.ndarray, eps: float | None = None):
        """(fwd_flux, rev_flux) per channel at concentration vector x."""
        if eps is None:
            eps = self.eps
        x = np.asarray(x, dtype=float)
        x_app = np.append(x, 1.0)
        Rf, Rr, _ = self._factors()
        kf, kr, gate_f, gate_r = self._rate_vectors()
        C = self.catalysis_matrix()
        cat = self.kappa * np.asarray(
            C @ (self.cat_weights
                 * np.power(x, self.cat_power))).ravel()
        if self.k_inhibit:
            # product inhibition / active-site saturation: the
            # channel's own catalyst pool suppresses marginal boost —
            # identical divisor on both directions, so k_f/k_r and
            # detailed balance are preserved
            cat = cat / (1.0 + self.k_inhibit
                         * np.asarray(self._catalysis_binary()
                                      @ x).ravel())
        boost_f = np.where(gate_f, eps, 1.0) + cat
        boost_r = np.where(gate_r, eps, 1.0) + cat
        mf = np.prod(x_app[Rf], axis=1)
        mr = np.prod(x_app[Rr], axis=1)
        return kf * mf * boost_f, kr * mr * boost_r

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
        kf, kr, gate_f, gate_r = self._rate_vectors()
        C = self.catalysis_matrix()
        bg_f = np.where(gate_f, eps, 1.0)
        bg_r = np.where(gate_r, eps, 1.0)
        Cx = np.asarray(
            C @ (self.cat_weights
                 * np.power(x, self.cat_power))).ravel()
        # inhibition denominator per channel: 1 + K * sum_{c in C_r} x_c
        # (catalyst concentrations — unweighted binary incidence)
        if self.k_inhibit:
            denom = (1.0 + self.k_inhibit
                     * np.asarray(self._catalysis_binary() @ x).ravel())
        else:
            denom = np.ones(self.n_channels)
        mf = np.prod(x_app[Rf], axis=1)
        mr = np.prod(x_app[Rr], axis=1)

        Df = np.zeros((nc, n))
        Db = np.zeros((nc, n))
        for r, ch in enumerate(self.channels):
            fac_f = ch.factors("fwd")
            fac_r = ch.factors("rev")
            cat_term = (self.kappa * Cx[r] / denom[r]
                        if self.k_inhibit else self.kappa * Cx[r])
            # d(monomial)/dx_i: sum over positions holding species i of the
            # product over the other positions.
            for pos, s in enumerate(fac_f):
                others = np.prod(
                    [x_app[fac_f[q]] for q in range(len(fac_f)) if q != pos],
                    axis=0,
                ) if len(fac_f) > 1 else 1.0
                Df[r, s] += kf[r] * (bg_f[r] + cat_term) * others
            for pos, s in enumerate(fac_r):
                others = np.prod(
                    [x_app[fac_r[q]] for q in range(len(fac_r)) if q != pos],
                    axis=0,
                ) if len(fac_r) > 1 else 1.0
                Db[r, s] += kr[r] * (bg_r[r] + cat_term) * others
            # catalysts contribute via boost: d(k*mass*kappa*x_c^p)/dx_c
            cat_idx = np.nonzero(C[r])[0]
            if len(cat_idx):
                # derivative of kappa * w_c * x_c^p / denom:
                # quotient rule — kappa*(w_c*p*x_c^(p-1)*denom -
                # Cx_r*K)/denom^2; at K=0 reduces to kappa*w_c*p*x_c^(p-1).
                # Raw power, matching the unclipped rate law in rates();
                # integer cat_power keeps this finite at any real x.
                edge_w = C[r, cat_idx]
                if hasattr(edge_w, "todense"):
                    edge_w = np.asarray(edge_w.todense())
                dpow = (self.cat_power * self.cat_weights[cat_idx]
                        * np.power(x[cat_idx], self.cat_power - 1)
                        * np.asarray(edge_w).ravel())
                if self.k_inhibit:
                    dpow = ((dpow * denom[r]
                             - Cx[r] * self.k_inhibit)
                            / denom[r] ** 2)
                Df[r, cat_idx] += kf[r] * mf[r] * self.kappa * dpow
                Db[r, cat_idx] += kr[r] * mr[r] * self.kappa * dpow
        return Df, Db

    def _cat_idx_padded(self) -> tuple[np.ndarray, np.ndarray, int]:
        """(nc x mc) padded catalyst index array (sentinel n_species)
        and parallel per-edge weight array (1.0 padding/defaults)."""
        if "catidx" not in self._cache:
            cats = [ch.catalysts for ch in self.channels]
            mc = max([len(c) for c in cats] + [1])
            R = np.full((self.n_channels, mc), self.n_species,
                        dtype=int)
            W = np.ones((self.n_channels, mc))
            for r, c in enumerate(cats):
                R[r, : len(c)] = c
                ch = self.channels[r]
                if ch.cat_w:
                    W[r, : len(c)] = ch.cat_w
            self._cache["catidx"] = (R, W, mc)
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
        kf, kr, gate_f, gate_r = self._rate_vectors()
        C = self.catalysis_matrix()
        bg_f = np.where(gate_f, eps, 1.0)
        bg_r = np.where(gate_r, eps, 1.0)
        Cx = np.asarray(
            C @ (self.cat_weights
                 * np.power(x, self.cat_power))).ravel()
        if self.k_inhibit:
            denom = (1.0 + self.k_inhibit
                     * np.asarray(self._catalysis_binary() @ x).ravel())
        else:
            denom = np.ones(self.n_channels)
        cat_term = self.kappa * Cx / denom
        mf = np.prod(x_app[Rf], axis=1)
        mr = np.prod(x_app[Rr], axis=1)

        rows_f, cols_f, data_f = [], [], []
        rows_b, cols_b, data_b = [], [], []
        rr = np.arange(nc)
        for pos in range(m):
            for R, xk, kvec, flux_m, bg, rows, cols, data in (
                    (Rf, x_app, kf, mf, bg_f, rows_f, cols_f, data_f),
                    (Rr, x_app, kr, mr, bg_r, rows_b, cols_b, data_b)):
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
                data.append((kvec * (bg + cat_term)
                             * others)[mask])
        Rc, Wc, mc = self._cat_idx_padded()
        for pos in range(mc):
            c = Rc[:, pos]
            mask = c < n
            if not mask.any():
                continue
            dpow = (self.cat_power
                    * self.cat_weights[np.where(mask, c, 0)]
                    * np.power(x[np.where(mask, c, 0)],
                               self.cat_power - 1)
                    * Wc[:, pos])
            if self.k_inhibit:
                # quotient rule on kappa*Cx/denom (denominator is
                # per-channel, shared across catalyst positions)
                dpow = ((dpow * denom - Cx * self.k_inhibit)
                        / denom ** 2)
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
        retention: float = 0.0,
        len_scale: float = 0.0,
    ) -> np.ndarray:
        """J = S (Df - Db) - diag(dilution_outflow).

        retention/len_scale define the CSTR outflow vector (see
        dilution_outflow); zeros reproduce the homogeneous-dilution
        Jacobian.
        """
        Df, Db = self.d_rates(x, eps)
        if retention and not (0.0 <= retention < 1.0):
            raise ValueError("retention must be in [0, 1)")
        outflow = self.dilution_outflow(dilution, retention, len_scale)
        if self.use_sparse:
            import scipy.sparse as sp
            return (self.stoichiometry() @ (Df - Db)
                    - sp.diags(outflow, format="csr"))
        return self.stoichiometry() @ (Df - Db) - np.diag(outflow)

    # ---------------- thermodynamics ----------------------------------
    def check_detailed_balance(self, tol: float = 1e-6) -> dict:
        """Verify k_f/k_r == exp(-g . nu) per channel; also verify that a
        closed system with these g values is at equilibrium (report the
        per-channel residual of K_eq consistency).

        Irreversible channels (k_f == 0 xor k_r == 0) cannot satisfy
        detailed balance on finite free energies. They are FAILED unless
        explicitly declared with Channel(driven=True) — a labeled
        driven-limit boundary process (e.g. benchmark source/sink)."""
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
            if ch.k_f <= 0 or ch.k_r <= 0:
                # an irreversible channel has infinite |dG|; it cannot
                # satisfy k_f/k_r = exp(-dG) on finite g
                rows.append(
                    {"channel": r, "dG": dG, "K_eq": K,
                     "k_f": ch.k_f, "k_r": ch.k_r,
                     "error": "irreversible_channel",
                     "declared_driven": bool(ch.driven),
                     "mechanism": ch.mechanism}
                )
                if not ch.driven:
                    ok = False
                continue
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
        internal (non-boundary) channels.

        Networks declared `params["open_system"]=True` model clamped
        reservoirs / source-sink channels that intentionally do not
        conserve letter mass (e.g. benchmarks.schlogl). Violations are
        still REPORTED in full, but `ok` passes — the declaration makes
        the deliberate open-system design explicit instead of a silent
        invariant break."""
        if self.comp is None:
            raise ValueError("no composition matrix")
        declared = bool(self.params.get("open_system", False))
        S = self.stoichiometry()
        res = (np.asarray(S.T @ self.comp.T).T
               if self.use_sparse else self.comp @ S)
        bad = np.argwhere(np.abs(res) > 0)
        return {
            "ok": bad.size == 0 or declared,
            "declared_open": declared,
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
            "cat_len_power": self.cat_len_power,
            "k_inhibit": self.k_inhibit,
            "inlet_species": sorted(self.inlet_species),
            "params": self.params,
            "channels": [
                {
                    "reactants": [list(t) for t in ch.reactants],
                    "products": [list(t) for t in ch.products],
                    "catalysts": list(ch.catalysts),
                    "cat_w": list(ch.cat_w),
                    "k_f": ch.k_f,
                    "k_r": ch.k_r,
                    "synth": ch.synth,
                    "mechanism": ch.mechanism,
                    "synth_dir": ch.synth_dir,
                    "driven": ch.driven,
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
                synth_dir=c.get("synth_dir", "both"),
                driven=c.get("driven", False),
                cat_w=tuple(c.get("cat_w") or ()),
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
            cat_len_power=d.get("cat_len_power", 0.0),
            k_inhibit=d.get("k_inhibit", 0.0),
            inlet_species=d.get("inlet_species"),
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
            cat_len_power=self.cat_len_power,
            k_inhibit=self.k_inhibit,
            inlet_species=self.inlet_species,
            params=self.params,
        )
        kw.update(overrides)
        return ReactionNetwork(self.species, list(channels), self.food, **kw)

    def without_catalysis(self) -> "ReactionNetwork":
        """Intervention: remove all catalysts (k_cat -> 0, K_eq unchanged)."""
        return self.with_channels(
            [
                Channel(ch.reactants, ch.products, (), ch.k_f, ch.k_r,
                        ch.synth, ch.mechanism, ch.synth_dir, ch.driven)
                for ch in self.channels
            ]
        )
