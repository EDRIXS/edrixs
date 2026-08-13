import numpy as np
import pytest
import edrixs
from edrixs.utils import (
    beta_to_kelvin, kelvin_to_beta, boltz_dist,
    UJ_to_UdJH, UdJH_to_UJ, UdJH_to_F0F2F4, UdJH_to_F0F2F4F6,
    F0F2F4_to_UdJH, F0F2F4_to_UJ, F0F2F4F6_to_UdJH,
    info_atomic_shell, case_to_shell_name, edge_to_shell_name,
    slater_integrals_name, rescale
)


# --- Temperature conversion ---

@pytest.mark.parametrize("T", [100.0, 300.0, 1000.0, 5000.0])
def test_beta_kelvin_roundtrip(T):
    """beta_to_kelvin and kelvin_to_beta are mutual inverses."""
    beta = kelvin_to_beta(T)
    T2 = beta_to_kelvin(beta)
    assert np.isclose(T, T2)


def test_kelvin_to_beta_scales_inversely():
    """Doubling temperature halves beta."""
    T = 300.0
    beta1 = kelvin_to_beta(T)
    beta2 = kelvin_to_beta(2 * T)
    assert np.isclose(beta1 / beta2, 2.0)


def test_beta_to_kelvin_known_value():
    """beta=1 corresponds to T = 1/(kB) ≈ 11604 K."""
    kb = 8.6173303e-5  # eV/K
    T = beta_to_kelvin(1.0)
    assert np.isclose(T, 1.0 / kb)


# --- Boltzmann distribution ---

@pytest.mark.parametrize("T", [100.0, 300.0, 1000.0])
def test_boltz_dist_normalized(T):
    """Boltzmann distribution sums to 1."""
    gs = [0.0, 0.1, 0.5, 1.0]
    dist = boltz_dist(gs, T)
    assert np.isclose(np.sum(dist), 1.0)


def test_boltz_dist_ground_state_dominant_at_low_T():
    """Ground state population approaches 1 at very low temperature."""
    gs = [0.0, 10.0]  # large gap ensures ground state dominance at low T
    dist = boltz_dist(gs, T=0.01)
    assert dist[0] > 0.999


def test_boltz_dist_equal_at_high_T():
    """At very high temperature, Boltzmann distribution becomes uniform."""
    gs = [0.0, 1.0, 2.0, 3.0]
    dist = boltz_dist(gs, T=1e8)
    expected = np.ones(4) / 4
    assert np.allclose(dist, expected, atol=1e-4)


def test_boltz_dist_ordering():
    """Lower energy levels have higher population."""
    gs = [0.0, 0.5, 1.0]
    dist = boltz_dist(gs, T=300.0)
    assert dist[0] > dist[1] > dist[2]


def test_boltz_dist_non_negative():
    """Boltzmann distribution values are all non-negative."""
    gs = [-1.0, 0.0, 1.0, 2.0]
    dist = boltz_dist(gs, T=300.0)
    assert np.all(dist >= 0)


# --- Kanamori/Slater integral conversions ---

@pytest.mark.parametrize("U,J", [(3.0, 0.5), (4.0, 0.8), (2.5, 0.3)])
def test_UJ_UdJH_roundtrip(U, J):
    """UJ_to_UdJH and UdJH_to_UJ are mutual inverses."""
    Ud, JH = UJ_to_UdJH(U, J)
    U2, J2 = UdJH_to_UJ(Ud, JH)
    assert np.isclose(U, U2)
    assert np.isclose(J, J2)


@pytest.mark.parametrize("Ud,JH", [(4.0, 0.8), (3.0, 0.5), (5.0, 1.0)])
def test_UdJH_F0F2F4_roundtrip(Ud, JH):
    """UdJH_to_F0F2F4 and F0F2F4_to_UdJH are mutual inverses."""
    F0, F2, F4 = UdJH_to_F0F2F4(Ud, JH)
    Ud2, JH2 = F0F2F4_to_UdJH(F0, F2, F4)
    assert np.isclose(Ud, Ud2)
    assert np.isclose(JH, JH2)


@pytest.mark.parametrize("Ud,JH", [(4.0, 0.8), (3.0, 0.5)])
def test_UdJH_F0F2F4F6_roundtrip(Ud, JH):
    """UdJH_to_F0F2F4F6 and F0F2F4F6_to_UdJH are mutual inverses."""
    F0, F2, F4, F6 = UdJH_to_F0F2F4F6(Ud, JH)
    Ud2, JH2 = F0F2F4F6_to_UdJH(F0, F2, F4, F6)
    assert np.isclose(Ud, Ud2)
    assert np.isclose(JH, JH2)


def test_UdJH_to_F0F2F4_F0_equals_Ud():
    """For d-orbitals, F0 = Ud."""
    Ud, JH = 4.0, 0.8
    F0, F2, F4 = UdJH_to_F0F2F4(Ud, JH)
    assert np.isclose(F0, Ud)


def test_UdJH_to_F0F2F4F6_F0_equals_Ud():
    """For f-orbitals, F0 = Ud."""
    Ud, JH = 5.0, 1.0
    F0, F2, F4, F6 = UdJH_to_F0F2F4F6(Ud, JH)
    assert np.isclose(F0, Ud)


def test_F0F2F4_to_UdJH_formula():
    """F0F2F4_to_UdJH: Ud=F0, JH=(F2+F4)/14."""
    F0, F2, F4 = 4.0, 10.0, 6.0
    Ud, JH = F0F2F4_to_UdJH(F0, F2, F4)
    assert np.isclose(Ud, F0)
    assert np.isclose(JH, (F2 + F4) / 14.0)


def test_F0F2F4_to_UJ_formula():
    """F0F2F4_to_UJ: U = F0 + 4/49*(F2+F4), J = 3/49*F2 + 20/441*F4."""
    F0, F2, F4 = 4.0, 10.0, 6.0
    U, J = F0F2F4_to_UJ(F0, F2, F4)
    expected_U = F0 + 4.0 / 49.0 * (F2 + F4)
    expected_J = 3.0 / 49.0 * F2 + 20 / 441.0 * F4
    assert np.isclose(U, expected_U)
    assert np.isclose(J, expected_J)


# --- info_atomic_shell ---

def test_info_atomic_shell_contains_expected_keys():
    """info_atomic_shell contains all expected shell names."""
    info = info_atomic_shell()
    for key in ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']:
        assert key in info


def test_info_atomic_shell_d_values():
    """d-shell has l=2 and 10 orbitals."""
    info = info_atomic_shell()
    assert info['d'][0] == 2
    assert info['d'][1] == 10


def test_info_atomic_shell_p_values():
    """p-shell has l=1 and 6 orbitals."""
    info = info_atomic_shell()
    assert info['p'][0] == 1
    assert info['p'][1] == 6


def test_info_atomic_shell_f_values():
    """f-shell has l=3 and 14 orbitals."""
    info = info_atomic_shell()
    assert info['f'][0] == 3
    assert info['f'][1] == 14


def test_info_atomic_shell_t2g_values():
    """t2g has l=2 (effective) and 6 orbitals."""
    info = info_atomic_shell()
    assert info['t2g'][0] == 2
    assert info['t2g'][1] == 6


# --- case_to_shell_name ---

@pytest.mark.parametrize("case,expected", [
    ('d', ('d',)),
    ('p', ('p',)),
    ('t2g', ('t2g',)),
    ('f', ('f',)),
    ('dp', ('d', 'p')),
    ('t2gp32', ('t2g', 'p32')),
    ('dp12', ('d', 'p12')),
    ('fd', ('f', 'd')),
])
def test_case_to_shell_name(case, expected):
    """case_to_shell_name returns correct shell name tuple."""
    assert case_to_shell_name(case) == expected


# --- edge_to_shell_name ---

@pytest.mark.parametrize("edge,expected_shell", [
    ('L3', 'p32'),
    ('L2', 'p12'),
    ('L23', 'p'),
    ('K', 's'),
    ('M4', 'd32'),
    ('M5', 'd52'),
    ('M45', 'd'),
    ('N6', 'f52'),
    ('N7', 'f72'),
    ('N67', 'f'),
])
def test_edge_to_shell_name(edge, expected_shell):
    """edge_to_shell_name returns correct shell name."""
    assert edge_to_shell_name(edge) == expected_shell


def test_edge_to_shell_name_with_main_qn():
    """edge_to_shell_name with_main_qn=True includes principal quantum number."""
    assert edge_to_shell_name('L3', with_main_qn=True) == '2p32'
    assert edge_to_shell_name('M45', with_main_qn=True) == '3d'
    assert edge_to_shell_name('N67', with_main_qn=True) == '4f'


# --- slater_integrals_name ---

def test_slater_integrals_name_d_single_shell():
    """d-shell has Slater integrals F0_11, F2_11, F4_11."""
    names = slater_integrals_name(('d',))
    assert names == ['F0_11', 'F2_11', 'F4_11']


def test_slater_integrals_name_p_single_shell():
    """p-shell has Slater integrals F0_11, F2_11."""
    names = slater_integrals_name(('p',))
    assert names == ['F0_11', 'F2_11']


def test_slater_integrals_name_f_single_shell():
    """f-shell has Slater integrals F0_11, F2_11, F4_11, F6_11."""
    names = slater_integrals_name(('f',))
    assert names == ['F0_11', 'F2_11', 'F4_11', 'F6_11']


def test_slater_integrals_name_dp_two_shells():
    """d+p two-shell case contains cross terms F2_12, G1_12, G3_12."""
    names = slater_integrals_name(('d', 'p'))
    assert 'F0_11' in names
    assert 'F2_11' in names
    assert 'F4_11' in names
    assert 'F0_12' in names
    assert 'F2_12' in names
    assert 'G1_12' in names
    assert 'G3_12' in names
    assert 'F0_22' in names
    assert 'F2_22' in names


def test_slater_integrals_name_with_custom_label():
    """Custom labels appear in Slater integral names."""
    names = slater_integrals_name(('d',), label=('d',))
    assert 'F0_dd' in names
    assert 'F2_dd' in names
    assert 'F4_dd' in names


def test_slater_integrals_name_dp_with_labels():
    """Custom labels for d+p shells produce labeled names."""
    names = slater_integrals_name(('d', 'p'), label=('d', 'p'))
    assert 'F0_dd' in names
    assert 'G1_dp' in names
    assert 'F2_pp' in names


# --- rescale ---

def test_rescale_with_no_scale():
    """rescale with scale=None returns unchanged list."""
    old_list = [1.0, 2.0, 3.0]
    new_list = rescale(old_list, scale=None)
    assert new_list == old_list


def test_rescale_applies_factors():
    """rescale multiplies specified positions by given factors."""
    old_list = [1.0, 2.0, 3.0, 4.0]
    new_list = rescale(old_list, ([1, 3], [2.0, 0.5]))
    assert np.isclose(new_list[0], 1.0)
    assert np.isclose(new_list[1], 4.0)   # 2.0 * 2.0
    assert np.isclose(new_list[2], 3.0)
    assert np.isclose(new_list[3], 2.0)   # 4.0 * 0.5


def test_rescale_does_not_modify_original():
    """rescale does not modify the original list."""
    old_list = [1.0, 2.0, 3.0]
    original_copy = old_list.copy()
    rescale(old_list, ([0], [10.0]))
    assert old_list == original_copy


def test_rescale_single_element():
    """rescale with a single scaling factor."""
    old_list = [5.0, 10.0, 15.0]
    new_list = rescale(old_list, ([0], [3.0]))
    assert np.isclose(new_list[0], 15.0)
    assert np.isclose(new_list[1], 10.0)
    assert np.isclose(new_list[2], 15.0)


# --- get_atom_data ---

def test_get_atom_data_ni_3d8_structure():
    """get_atom_data for Ni 3d8 returns expected Slater integral structure."""
    res = edrixs.get_atom_data('Ni', v_name='3d', v_noccu=8)
    assert 'slater_i' in res
    assert 'v_soc_i' in res
    names = [x[0] for x in res['slater_i']]
    assert 'F0_11' in names
    assert 'F2_11' in names
    assert 'F4_11' in names


def test_get_atom_data_ni_with_edge():
    """get_atom_data with edge returns both initial and intermediate Hamiltonians."""
    res = edrixs.get_atom_data('Ni', v_name='3d', v_noccu=8, edge='L3')
    assert 'slater_i' in res
    assert 'slater_n' in res
    assert 'v_soc_i' in res
    assert 'v_soc_n' in res
    assert 'c_soc' in res
    assert 'edge_ene' in res
    assert 'gamma_c' in res


def test_get_atom_data_soc_nonzero():
    """SOC for Ni 3d shell is non-zero."""
    res = edrixs.get_atom_data('Ni', v_name='3d', v_noccu=8)
    assert res['v_soc_i'][0] > 0


# --- Independent regression oracles added after review ---

def test_boltz_dist_two_level_exact_ratio():
    """Pin the quantitative Boltzmann factor, not only normalization/order."""
    energies = [0.0, 0.25]
    temperature = 300.0
    kb = 8.6173303e-5
    dist = boltz_dist(energies, temperature)
    assert np.isclose(dist[1] / dist[0], np.exp(-0.25 / (kb * temperature)))


def test_UJ_to_UdJH_matches_explicit_formula():
    """An explicit formula prevents mutually-wrong inverse functions from passing."""
    U, J = 4.0, 0.8
    F2 = J / (3.0 / 49.0 + 20.0 * 0.625 / 441.0)
    F4 = 0.625 * F2
    expected = (U - 4.0 / 49.0 * (F2 + F4), (F2 + F4) / 14.0)
    assert np.allclose(UJ_to_UdJH(U, J), expected)


def test_UdJH_to_F0F2F4_pins_F4_over_F2_ratio():
    """Roundtrip tests alone do not constrain the conventional F4/F2 ratio."""
    F0, F2, F4 = UdJH_to_F0F2F4(4.0, 0.8)
    assert np.isclose(F0, 4.0)
    assert np.isclose(F4 / F2, 0.625)


def test_UdJH_to_F0F2F4F6_pins_standard_ratios():
    """Pin both f-shell Slater-integral ratios used by the conversion."""
    F0, F2, F4, F6 = UdJH_to_F0F2F4F6(5.0, 1.0)
    assert np.isclose(F0, 5.0)
    assert np.isclose(F4 / F2, 451.0 / 675.0)
    assert np.isclose(F6 / F2, 1001.0 / 2025.0)
