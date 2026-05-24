import numpy as np
import pytest
from math import comb
import edrixs
from edrixs.fock_basis import (
    fock_bin, get_fock_bin_by_N, get_fock_full_N,
    get_fock_basis_by_NLz, get_fock_basis_by_NSz,
    get_fock_basis_by_NJz, get_fock_basis_by_N_LzSz
)


@pytest.mark.parametrize("n,k", [(4, 2), (5, 3), (6, 0), (6, 6), (3, 1), (5, 0)])
def test_fock_bin_count(n, k):
    """fock_bin(n, k) produces C(n, k) states."""
    result = fock_bin(n, k)
    assert len(result) == comb(n, k)


@pytest.mark.parametrize("n,k", [(4, 2), (5, 3), (3, 1), (4, 4)])
def test_fock_bin_each_state_has_correct_length_and_occupancy(n, k):
    """Each state from fock_bin has length n and exactly k ones."""
    for state in fock_bin(n, k):
        assert len(state) == n
        assert sum(state) == k


def test_fock_bin_n4_k2_exact():
    """fock_bin(4, 2) matches the canonical set of states."""
    result = fock_bin(4, 2)
    expected = [
        [1, 1, 0, 0],
        [1, 0, 1, 0],
        [1, 0, 0, 1],
        [0, 1, 1, 0],
        [0, 1, 0, 1],
        [0, 0, 1, 1],
    ]
    assert result == expected


def test_fock_bin_n4_k0_is_empty_state():
    """fock_bin(4, 0) gives a single all-zero state."""
    result = fock_bin(4, 0)
    assert result == [[0, 0, 0, 0]]


def test_get_fock_bin_by_N_single_shell_matches_fock_bin():
    """get_fock_bin_by_N with one shell is identical to fock_bin."""
    assert get_fock_bin_by_N(4, 2) == fock_bin(4, 2)
    assert get_fock_bin_by_N(5, 3) == fock_bin(5, 3)


def test_get_fock_bin_by_N_two_shells_count():
    """Two-shell basis has C(n1,k1) × C(n2,k2) states."""
    n1, k1, n2, k2 = 4, 2, 2, 1
    result = get_fock_bin_by_N(n1, k1, n2, k2)
    assert len(result) == comb(n1, k1) * comb(n2, k2)


def test_get_fock_bin_by_N_two_shells_length_and_occupancy():
    """Each two-shell state has correct total length and total occupancy."""
    n1, k1, n2, k2 = 4, 2, 2, 1
    for state in get_fock_bin_by_N(n1, k1, n2, k2):
        assert len(state) == n1 + n2
        assert sum(state) == k1 + k2


def test_get_fock_full_N_known_values():
    """get_fock_full_N(4, 2) returns the 6 expected decimal Fock states."""
    result = sorted(get_fock_full_N(4, 2))
    # Binary: 0011=3, 0101=5, 0110=6, 1001=9, 1010=10, 1100=12
    expected = sorted([3, 5, 6, 9, 10, 12])
    assert result == expected


@pytest.mark.parametrize("norb,N", [(4, 0), (4, 2), (4, 4), (6, 3), (8, 4)])
def test_get_fock_full_N_count(norb, N):
    """get_fock_full_N returns C(norb, N) states."""
    result = get_fock_full_N(norb, N)
    assert len(result) == comb(norb, N)


def test_get_fock_full_N_each_state_has_N_bits():
    """Each state from get_fock_full_N has exactly N bits set."""
    norb, N = 6, 3
    for state in get_fock_full_N(norb, N):
        assert bin(state).count('1') == N


def test_get_fock_basis_by_NLz_total_count():
    """Total number of states across all Lz sectors equals C(norb, N)."""
    norb, N = 6, 2
    lz_list = [-1, -1, 0, 0, 1, 1]
    basis = get_fock_basis_by_NLz(norb, N, lz_list)
    total = sum(len(v) for v in basis.values())
    assert total == comb(norb, N)


def test_get_fock_basis_by_NLz_extreme_sectors():
    """Lz = ±(max_lz * N) each contain exactly one state."""
    norb, N = 6, 2
    lz_list = [-1, -1, 0, 0, 1, 1]
    basis = get_fock_basis_by_NLz(norb, N, lz_list)
    # Lz = -2: only state where both electrons are in m=-1 orbitals
    assert len(basis[-2]) == 1
    # Lz = +2: only state where both electrons are in m=+1 orbitals
    assert len(basis[2]) == 1


def test_get_fock_basis_by_NSz_total_count():
    """Total number of states across all Sz sectors equals C(norb, N)."""
    norb, N = 6, 2
    sz_list = [1, -1, 1, -1, 1, -1]
    basis = get_fock_basis_by_NSz(norb, N, sz_list)
    total = sum(len(v) for v in basis.values())
    assert total == comb(norb, N)


def test_get_fock_basis_by_NSz_keys_span_correct_range():
    """Sz basis has keys from min_sz*N to max_sz*N."""
    norb, N = 4, 2
    sz_list = [1, -1, 1, -1]
    basis = get_fock_basis_by_NSz(norb, N, sz_list)
    assert -2 in basis  # Sz = -1 - 1 = -2
    assert 2 in basis   # Sz = +1 + 1 = +2
    assert 0 in basis


def test_get_fock_basis_by_NJz_total_count():
    """Total number of states across all Jz sectors equals C(norb, N)."""
    norb, N = 6, 2
    jz_list = [-1, 1, -3, -1, 1, 3]
    basis = get_fock_basis_by_NJz(norb, N, jz_list)
    total = sum(len(v) for v in basis.values())
    assert total == comb(norb, N)


def test_get_fock_basis_by_N_LzSz_total_count():
    """Total number of states across all (Lz, Sz) sectors equals C(norb, N)."""
    norb, N = 6, 2
    lz_list = [-1, -1, 0, 0, 1, 1]
    sz_list = [1, -1, 1, -1, 1, -1]
    basis = get_fock_basis_by_N_LzSz(norb, N, lz_list, sz_list)
    total = sum(len(v) for v in basis.values())
    assert total == comb(norb, N)


def test_get_fock_basis_by_N_LzSz_keys_are_tuples():
    """Keys in (Lz, Sz) basis are tuples of two integers."""
    norb, N = 4, 1
    lz_list = [-1, -1, 0, 0]
    sz_list = [1, -1, 1, -1]
    basis = get_fock_basis_by_N_LzSz(norb, N, lz_list, sz_list)
    for key in basis.keys():
        assert isinstance(key, tuple)
        assert len(key) == 2


def test_get_fock_basis_by_NLz_states_have_correct_Lz():
    """Each state in a given Lz sector actually has that Lz value."""
    norb, N = 6, 2
    lz_list = [-1, -1, 0, 0, 1, 1]
    basis = get_fock_basis_by_NLz(norb, N, lz_list)
    for Lz_key, states in basis.items():
        for state in states:
            # Compute Lz from the decimal representation
            Lz_computed = sum(lz_list[i] for i in range(norb) if (state >> i) & 1)
            assert Lz_computed == Lz_key
