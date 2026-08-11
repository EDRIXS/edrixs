"""Unit tests for SciPy sparse many-body operator construction."""

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

from edrixs.fock_basis import FockBasis
from edrixs.scipy_backend.scipy_backend import (
    _four_fermion_csr_from_sparse_umat,
    four_fermion_csr,
    four_fermion_csr_auto,
    two_fermion_csr,
)
from edrixs.solvers import build_op

from _oracles import fixed_particle_basis, one_body_oracle, two_body_oracle


def test_fock_basis_encode_decode_roundtrip():
    """Fock-state encodings round-trip through their basis positions."""
    basis = FockBasis([0b1100, 0b1010, 0b1001], norbs=4)

    assert len(basis) == 3
    for position, state in enumerate(basis.basis_int):
        assert basis.encode(state) == position
        assert basis.decode(position) == state


def test_one_body_csr_matches_independent_jordan_wigner_oracle():
    """The SciPy one-body constructor matches an independent oracle."""
    rng = np.random.default_rng(4)
    basis = fixed_particle_basis(4, 2)
    emat = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))

    actual = two_fermion_csr(emat, basis).toarray()
    expected = one_body_oracle(emat, basis)

    assert_allclose(actual, expected, rtol=0, atol=1e-13)


def test_one_particle_sector_reproduces_orbital_matrix():
    """A one-particle Fock-space operator reproduces its orbital matrix."""
    basis = FockBasis([0b10, 0b01], norbs=2)
    emat = np.array([[1.2, 0.3 + 0.4j], [0.3 - 0.4j, -0.2]])

    assert_allclose(two_fermion_csr(emat, basis).toarray(), emat)


def test_four_body_dense_and_flat_sparse_paths_match_oracle():
    """Dense and flattened sparse Coulomb inputs match an independent oracle."""
    basis = fixed_particle_basis(4, 2)
    umat = np.zeros((4, 4, 4, 4), dtype=complex)
    umat[0, 1, 1, 0] = 1.7
    umat[3, 2, 1, 0] = -0.2 + 0.3j
    umat[0, 1, 2, 3] = 0.4 - 0.1j

    expected = two_body_oracle(umat, basis)
    dense_path = four_fermion_csr(umat, basis).toarray()
    flat_sparse_u = sp.csr_matrix(umat.reshape(16, 16))
    sparse_path = four_fermion_csr_auto(flat_sparse_u, basis).toarray()

    assert_allclose(dense_path, expected, rtol=0, atol=1e-13)
    assert_allclose(sparse_path, expected, rtol=0, atol=1e-13)


def test_public_build_op_matches_backend_operator_construction():
    """The public ``build_op`` wrapper delegates to the SciPy constructor."""
    basis = fixed_particle_basis(3, 2)
    emat = np.diag([0.2, 0.5, 1.1])
    umat = np.zeros((3, 3, 3, 3), dtype=complex)
    umat[0, 1, 1, 0] = 0.7

    actual = build_op(emat, umat, basis, backend="scipy")
    expected = two_fermion_csr(emat, basis) + four_fermion_csr(umat, basis)

    assert_allclose(actual.toarray(), expected.toarray())


def test_tolerance_drops_small_terms():
    """The sparse constructor prunes only terms below the selected tolerance."""
    basis = fixed_particle_basis(2, 1)
    emat = np.diag([1e-11, 2e-9])

    actual = two_fermion_csr(emat, basis, tol=1e-10).toarray()

    assert actual[0, 0] == 0
    assert actual[1, 1] == pytest.approx(2e-9)


def test_basis_norbs_mismatch_raises():
    """Transition construction rejects bases from different orbital spaces."""
    left = FockBasis([0b10, 0b01], norbs=2)
    right = FockBasis([0b100], norbs=3)

    with pytest.raises(ValueError, match="same norbs"):
        two_fermion_csr(np.eye(2), left, right)


def test_sparse_u_shape_is_validated():
    """Flattened sparse interactions require one axis per orbital pair."""
    basis = fixed_particle_basis(3, 2)

    with pytest.raises(ValueError, match="expected"):
        _four_fermion_csr_from_sparse_umat(sp.eye(8), basis)
