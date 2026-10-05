"""Analytic benchmark networks (test fixtures only — never headline models).

transcritical_motif
    The plan-§2.1 minimal autocatalytic motif: A + B <=> C catalyzed by C
    in a CSTR with inlet concentration F of A and B and dilution D
    (feed term D*(inlet - x); dead state has a = b = F). With eps -> 0
    and k_r -> 0 the steady-state condition factorizes and there is an
    exact transcritical exchange of stability at

        lambda_c :  k_f * kappa * F^2 = D   i.e.  k_f*kappa = D/F^2

    below lambda_c the dead state C=0 is stable and the other positive
    root is unphysical; above it C=0 is unstable and one physical
    organized root appears, growing as O ~ (lambda - lambda_c) (beta = 1).
    At eps > 0 the transition rounds into the imperfect-bifurcation form

        O = O0 + (A/2) [ (lambda - lambda_c)
                         + sqrt((lambda - lambda_c)^2 + 4 c eps) ].

schlogl
    Reversible driven Schlogl-type network: A + 2X <=> 3X and B <=> X in a
    CSTR. Provides a genuine fold/bistable regime for validating
    hysteresis/continuation machinery.
"""

from __future__ import annotations

import numpy as np

from .core import Channel, ReactionNetwork


def transcritical_motif(
    f_inlet: float = 1.0,
    k_lig: float = 1.0,
    kappa: float = 1.0,
    keq: float = np.inf,
    eps: float = 0.0,
) -> tuple[ReactionNetwork, dict]:
    """Return (network, analytic) for the A + B <=> C motif.

    analytic = {"lambda_c": k_f*kappa at the transcritical point,
                "c_star": function(lambda_kfkappa, D) -> physical root or 0}.
    """
    species = ("A", "B", "C")
    food = (0, 1)
    food_inlet = np.array([f_inlet, f_inlet, 0.0])
    g = np.array([0.0, 0.0, -np.log(keq) if np.isfinite(keq) else 0.0])
    comp = np.array([[1.0, 1.0, 2.0]])
    k_r = k_lig / keq if np.isfinite(keq) and keq > 0 else 0.0
    # keq = inf (irreversible) is a labeled edge case: use a finite
    # effective keq via g so reversibility machinery still applies.
    if not np.isfinite(keq):
        # k_r -> 0 limit; keep channel two-way-capable by storing a tiny
        # but nonzero effective ratio? No: k_r = 0 is a legitimate
        # one-way driven limit, flagged in params.
        pass
    ch = Channel(
        reactants=((0, 1), (1, 1)),
        products=((2, 1),),
        catalysts=(2,),
        k_f=k_lig,
        k_r=k_r,
        synth=True,
        mechanism="elongation",
    )
    net = ReactionNetwork(
        species, [ch], food,
        structure={2: "C"},
        comp=comp,
        g=g,
        food_inlet=food_inlet,
        kappa=kappa,
        eps=eps,
        params={"generator": "transcritical_motif", "f_inlet": f_inlet,
                "k_lig": k_lig, "kappa": kappa, "keq": keq, "eps": eps},
    )

    def c_star(lam: float, D: float, K_eq: float = keq) -> float:
        """Physical organized root at eps=0, or 0 if dead state is the
        only physical equilibrium. lam = k_f*kappa (the control); the
        fixed-point equation is k_f (F-c)^2 - k_r c = D/kappa, i.e.

            k_f c^2 - (2 k_f F + k_r) c + k_f F^2 - D/kappa = 0.
        """
        F, D_ = f_inlet, D
        k_r_ = k_lig / K_eq if np.isfinite(K_eq) and K_eq > 0 else 0.0
        kappa_ = lam / k_lig
        A = k_lig
        B = -(2 * k_lig * F + k_r_)
        C = k_lig * F * F - D_ / kappa_
        disc = B * B - 4 * A * C
        if disc < 0:
            return 0.0
        r = (-B - np.sqrt(disc)) / (2 * A)
        return float(r) if 0 < r < F else 0.0

    # lambda_c depends on D (a reactor parameter), so expose functions.
    analytic = {
        "lambda_c_of_D": lambda D: D / f_inlet ** 2,
        "kappa_c_of_D": lambda D: D / (k_lig * f_inlet ** 2),
        "c_star": c_star,
    }
    return net, analytic


def schlogl(
    a_in: float = 0.5,
    b_in: float = 0.5,
    k1: float = 6.0,
    k2: float = 1.0,
    k3: float = 6.0,
    k4: float = 65.0 / 6.0,
    kappa: float = 10.0,
) -> ReactionNetwork:
    """Driven Schlogl-type network with clamped (catalyst-role) reservoirs.

    A and B are food species that modulate rates catalytically but are
    not consumed by the channels, so their steady-state concentrations
    are simply their inlets. Both directions of each channel share the
    boost beta = 1 + kappa*x_reservoir = 6, so with the defaults

        dx/dt = beta*(k1 x^2 - k2 x^3) + beta*(k3 - k4 x) - D x
              = 6*(6 x^2 - x^3) + 6*(6 - 65/6 x) - D x
              = -6 (x-1)(x-2)(x-3) + (1 - D) x

    which at D=1 is exactly -6(x-1)(x-2)(x-3): stable roots at x=1 and
    x=3, unstable at x=2. D sweeps trace a fold. Tests verify
    bistability numerically, not exact roots.
    """
    species = ("B", "A", "X")
    food = (0, 1)
    food_inlet = np.array([b_in, a_in, 0.0])
    comp = np.array([[1.0, 1.0, 1.0]])
    g = np.zeros(3)
    channels = [
        # 2X <=> 3X, both directions boosted by reservoir A
        Channel(
            reactants=((2, 2),),
            products=((2, 3),),
            catalysts=(1,),
            k_f=k1, k_r=k2,
            synth=False, mechanism="autocat",
        ),
        # source/sink <=> X boosted by reservoir B
        Channel(
            reactants=(),
            products=((2, 1),),
            catalysts=(0,),
            k_f=k3, k_r=k4,
            synth=False, mechanism="reservoir",
        ),
    ]
    return ReactionNetwork(
        species, channels, food,
        comp=comp, g=g, food_inlet=food_inlet, kappa=kappa, eps=1.0,
        params={"generator": "schlogl", "a_in": a_in, "b_in": b_in,
                "k1": k1, "k2": k2, "k3": k3, "k4": k4,
                "kappa": kappa},
    )


def immigration_death(k_birth: float = 2.0, k_death: float = 0.5,
                      f_inlet: float = 0.0) -> ReactionNetwork:
    """X birth (from buffered food F -> X) + death X -> F.

    Stationary distribution of n_X in a closed... actually open system:
    birth is Poisson immigration (rate k_birth * x_F) and death is
    k_death * n_X — a classic testbed for SSA statistics.
    """
    species = ("F", "X")
    food = (0,)
    food_inlet = np.array([f_inlet, 0.0])
    comp = np.array([[1.0, 1.0]])
    g = np.zeros(2)
    channels = [
        Channel(reactants=((0, 1),), products=((1, 1),),
                k_f=k_birth, k_r=k_death, synth=False,
                mechanism="immigration_death"),
    ]
    return ReactionNetwork(
        species, channels, food, comp=comp, g=g,
        food_inlet=food_inlet, eps=1.0,
        params={"generator": "immigration_death", "k_birth": k_birth,
                "k_death": k_death},
    )
