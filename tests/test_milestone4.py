import numpy as np

from chemtrans.dynamics.ode import simulate
from chemtrans.measure import measure_run
from chemtrans.networks.generators import polymer_ligation
from chemtrans.observables.assembly import (
    _string_assembly,
    molecular_assembly_index,
)
from chemtrans.observables.heredity import (
    heredity_index,
    template_fidelity,
)


def test_string_assembly_index():
    m = {}
    assert _string_assembly("A", m) == 0
    assert _string_assembly("AA", m) == 1
    assert _string_assembly("AB", m) == 1
    assert _string_assembly("AAAA", m) == 2      # AA + AA, reuse subassembly
    assert _string_assembly("ABAB", m) == 2      # AB + AB
    assert _string_assembly("ABC", m) == 2       # (A+BC) or (AB+C)
    # repeats are cheaper than distinct strings of same length
    assert _string_assembly("AAAAAAAA", m) < _string_assembly("ABCDEFGH", m)


def test_polymer_generator_structure():
    net = polymer_ligation(alphabet_size=3, max_length=8, n_ligations=60,
                           p_cat=0.5, p_info=1.0, seed=4)
    assert len(net.food) == 3                       # monomers
    assert all(len(net.structure[i]) == 1 for i in net.food)
    # species are strings; longer ones exist
    lens = [len(s) for s in net.structure.values()]
    assert max(lens) > 1
    # reactions are ligation (2 reactants) or cleavage (1 reactant)
    assert all(len(r.reactants) <= 2 for r in net.reactions)
    # template catalysts contain their product
    templ = 0
    for rxn in net.reactions:
        for c in rxn.catalysts:
            for j in rxn.products:
                assert (net.structure[j] in net.structure[c]
                        or net.structure[c] in net.structure[j])
                templ += 1
    assert templ > 0


def test_molecular_assembly_index():
    net = polymer_ligation(alphabet_size=2, max_length=6, n_ligations=40,
                           seed=2)
    A = molecular_assembly_index(net)
    assert A[0] == 0 and A[1] == 0                # monomers
    assert np.isfinite(A).all()                   # every string has an index
    assert A.max() > 0


def test_template_fidelity_tracks_p_info():
    net_hi = polymer_ligation(alphabet_size=3, max_length=7, n_ligations=50,
                              p_cat=0.6, p_info=1.0, seed=1)
    net_lo = polymer_ligation(alphabet_size=3, max_length=7, n_ligations=50,
                              p_cat=0.6, p_info=0.0, seed=1)
    assert template_fidelity(net_hi) == 1.0
    assert template_fidelity(net_lo) < 0.5


def test_polymer_simulation_and_measure():
    net = polymer_ligation(alphabet_size=3, max_length=7, n_ligations=60,
                           p_cat=0.4, k_sat=1.0, seed=9)
    row = measure_run(net, feed_rate=1.0, dilution_rate=0.1, t_max=80)
    assert "H" in row
    assert 0.0 <= row["H"] <= 1.0
    assert row["assembly"]["max"] >= 1
