"""Generator invariants: mass conservation, thermodynamic consistency,
non-food catalysts, absorbing limit, deterministic hashing."""

import numpy as np

from chemtrans.networks.generators import polymer_network


def test_letter_conservation_all_channels():
    net = polymer_network(alphabet_size=2, max_len=3, seed=0)
    rep = net.check_letter_conservation()
    assert rep["ok"], rep["violations"]


def test_detailed_balance_by_construction():
    for dm in (0.0, 2.0, 5.0):
        net = polymer_network(alphabet_size=2, max_len=3, seed=1,
                              delta_mu=dm)
        assert net.check_detailed_balance()["ok"], dm


def test_catalysts_are_nonfood_polymers():
    net = polymer_network(alphabet_size=2, max_len=3, seed=2, p_cat=1.0,
                          min_cat_len=2)
    for ch in net.channels:
        for c in ch.catalysts:
            assert c not in net.food
            assert "*" not in net.species[c]
            assert len(net.species[c]) >= 2


def test_p_cat_zero_means_no_catalysis():
    net = polymer_network(alphabet_size=2, max_len=3, seed=3, p_cat=0.0)
    assert all(not ch.catalysts for ch in net.channels)


def test_p_selfrep_assigns_own_product():
    net = polymer_network(alphabet_size=2, max_len=3, seed=4, p_cat=1.0,
                          p_selfrep=1.0, min_cat_len=2)
    selfrep = [ch for ch in net.channels
               if ch.mechanism == "elongation" and ch.catalysts]
    assert selfrep
    for ch in selfrep:
        prod = {i for i, _ in ch.products}
        assert set(ch.catalysts) <= prod


def test_feed_mass_conserving_across_chi():
    a = polymer_network(alphabet_size=2, max_len=2, chi=0.0, seed=5)
    b = polymer_network(alphabet_size=2, max_len=2, chi=1.0, seed=5)
    assert np.isclose(a.feed_letters(), b.feed_letters())


def test_eps_zero_absorbing():
    """eps=0 + non-food catalysts => the zero-polymer state is absorbing
    for polymer synthesis (decay channels may still run)."""
    net = polymer_network(alphabet_size=2, max_len=3, seed=6, eps=0.0,
                          p_cat=0.5)
    x = net.food_inlet.copy()
    f, b = net.rates(x, eps=0.0)
    for r, ch in enumerate(net.channels):
        # each GATED direction must vanish at eps=0; an ungated decay
        # direction (thermal hydrolysis, synth_dir="rev" -> fwd ungated)
        # legitimately stays off here because its polymer reactant is
        # absent at the food-only state
        if ch.synth and ch.synth_dir != "rev":
            assert f[r] == 0.0
        if ch.synth and ch.synth_dir != "fwd":
            assert b[r] == 0.0


def test_eps_zero_absorbing_with_polymers_present():
    """The real absorbing-limit invariant: at eps=0, UNcatalyzed
    synthesis is zero even when polymer reactants are present — the
    food-only test above cannot catch a broken eps gate because synth
    monomials are zero there anyway."""
    net = polymer_network(alphabet_size=2, max_len=3, seed=6, eps=0.0,
                          p_cat=0.5)
    nonfood = [i for i in range(net.n_species) if not net.is_food[i]]
    # all polymers present, but every species that acts as a catalyst
    # for a synth channel is zeroed -> boost = eps + kappa*0 = 0
    x_nc = net.food_inlet.copy()
    x_nc[nonfood] = 0.5
    for ch in net.channels:
        for c in ch.catalysts:
            x_nc[c] = 0.0
    f, b = net.rates(x_nc, eps=0.0)
    for r, ch in enumerate(net.channels):
        if ch.synth and ch.synth_dir != "rev":
            assert f[r] == 0.0, \
                f"gated fwd of channel {r} fired at eps=0"
        if ch.synth and ch.synth_dir != "fwd":
            assert b[r] == 0.0, \
                f"gated rev of channel {r} fired at eps=0"
    # positive control: the same state at eps=1 must have nonzero
    # uncatalyzed synthesis (reactants are populated)
    f1, b1 = net.rates(x_nc, eps=1.0)
    assert any(ch.synth and (f1[r] > 0 or b1[r] > 0)
               for r, ch in enumerate(net.channels))
    # and at eps=0 a catalyzed synth channel CAN fire when its catalyst
    # is present
    net_cat = polymer_network(alphabet_size=2, max_len=3, seed=6,
                              eps=0.0, p_cat=1.0)
    cat_syn = [r for r, ch in enumerate(net_cat.channels)
               if ch.synth and ch.catalysts]
    assert cat_syn, "expected at least one catalyzed synth channel"
    r0 = cat_syn[0]
    x_c = net_cat.food_inlet.copy()
    for i, _ in net_cat.channels[r0].reactants:
        x_c[i] = 0.5
    for i, _ in net_cat.channels[r0].products:
        x_c[i] = 0.5
    x_c[net_cat.channels[r0].catalysts[0]] = 0.5
    fc, _ = net_cat.rates(x_c, eps=0.0)
    assert fc[r0] > 0.0


def test_thermal_hydrolysis_runs_at_eps_zero():
    """Invariant 3: decay channels still run at eps=0. Thermal
    hydrolysis (PA -> P + A, fwd) must fire whenever a polymer is
    present, while condensation (P + A -> PA, rev) stays eps-gated —
    synthesis needs the drive, cleavage does not. chi=0 feeds bare
    monomers so the condensation monomial is nonzero and the rev gate
    is genuinely exercised."""
    net = polymer_network(alphabet_size=2, max_len=3, seed=6, eps=0.0,
                          p_cat=0.5, chi=0.0)
    thermal = [r for r, ch in enumerate(net.channels)
               if ch.mechanism == "thermal"]
    assert thermal, "expected thermal channels"
    nonfood = [i for i in range(net.n_species) if not net.is_food[i]]
    x = net.food_inlet.copy()
    x[nonfood] = 0.5
    f, b = net.rates(x, eps=0.0)
    fired = sum(1 for r in thermal if f[r] > 0.0)
    assert fired > 0, "thermal hydrolysis frozen at eps=0"
    assert all(b[r] == 0.0 for r in thermal), \
        "thermal condensation escaped the eps=0 gate"
    # dead state stays absorbing: food-only state has no polymer for
    # hydrolysis to cleave and condensation is gated
    f0, b0 = net.rates(net.food_inlet.copy(), eps=0.0)
    assert all(f0[r] == 0.0 and b0[r] == 0.0 for r in thermal)
    # SSA agrees: hydrolysis processes exist with positive rate at eps=0
    from chemtrans.dynamics import gillespie
    procs = gillespie.build_processes(net, dilution=0.0, eps=0.0)
    th_fwd = [p for p in procs
              if p.kind == "channel" and p.channel in thermal
              and p.rate > 0]
    assert th_fwd, "SSA has no active thermal hydrolysis at eps=0"


def test_deterministic_hash():
    a = polymer_network(alphabet_size=2, max_len=3, seed=7)
    b = polymer_network(alphabet_size=2, max_len=3, seed=7)
    c = polymer_network(alphabet_size=2, max_len=3, seed=8)
    assert a.hash == b.hash
    assert a.hash != c.hash


def test_roundtrip_serialization():
    from chemtrans.networks.core import ReactionNetwork
    net = polymer_network(alphabet_size=2, max_len=3, seed=9)
    back = ReactionNetwork.from_dict(net.to_dict())
    assert back.hash == net.hash


def test_selfrep_respects_min_cat_len():
    """A product shorter than min_cat_len may not catalyze its own
    formation — no self-rep length leak."""
    net = polymer_network(alphabet_size=2, max_len=4, seed=0,
                          p_cat=1.0, p_selfrep=1.0, min_cat_len=3)
    assert any(ch.catalysts for ch in net.channels)
    for ch in net.channels:
        for c in ch.catalysts:
            assert len(net.species[c]) >= 3


def test_catalytic_multiplicity():
    """n_cat_max/p_multi add extra catalyst slots; defaults preserve the
    single-catalyst channel list exactly."""
    net = polymer_network(alphabet_size=2, max_len=3, seed=5,
                          p_cat=1.0, n_cat_max=3, p_multi=1.0)
    assert any(len(ch.catalysts) > 1 for ch in net.channels)
    for ch in net.channels:
        assert len(set(ch.catalysts)) == len(ch.catalysts)
    a = polymer_network(alphabet_size=2, max_len=3, seed=6)
    b = polymer_network(alphabet_size=2, max_len=3, seed=6,
                        n_cat_max=1, p_multi=0.0)
    assert [ch.catalysts for ch in a.channels] == \
           [ch.catalysts for ch in b.channels]


def test_ligation_channels():
    """P1+P2<=>P1P2: conservation, balance, no W, len>=2 substrates,
    absorbing limit, ordered-pair dedup."""
    net = polymer_network(alphabet_size=2, max_len=5, ligation=True,
                          p_cat=0.5, seed=1)
    lig = [ch for ch in net.channels if ch.mechanism == "ligation"]
    assert lig                       # channels exist at L=5
    iW = net.species.index("W")
    for ch in lig:
        assert all(i != iW for i, _ in ch.reactants + ch.products)
        assert all(len(net.species[i]) >= 2 for i, _ in ch.reactants)
        assert all(len(net.species[i])
                   == sum(a * len(net.species[j])
                          for j, a in ch.reactants)
                   for i, a in ch.products)
    assert net.check_letter_conservation()["ok"]
    assert net.check_detailed_balance()["ok"]
    # no duplicate (reactants, products, catalysts) signatures
    sigs = [(ch.reactants, ch.products, ch.catalysts) for ch in lig]
    assert len(sigs) == len(set(sigs))
    # eps=0 absorbing for ligation too (needs two existing polymers)
    net0 = polymer_network(alphabet_size=2, max_len=4, ligation=True,
                           eps=0.0, p_cat=0.5, seed=1)
    x = net0.food_inlet.copy()
    f, b = net0.rates(x, eps=0.0)
    for r, ch in enumerate(net0.channels):
        if ch.synth:
            assert f[r] == 0.0 and b[r] == 0.0
