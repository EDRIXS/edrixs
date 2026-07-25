"""Unit tests for sparse many-body operator construction.

These tests cover the stage between ``setup_*`` and the SciPy Krylov solvers.
They check that orbital-space basis data, one-body terms, and Coulomb terms are
converted into correct sparse Fock-space operators before ED, XAS, or RIXS.
"""

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

from edrixs.manybody_operator_csr import (
    FockBasis,
    _four_fermion_csr_from_sparse_umat,
    four_fermion_csr,
    four_fermion_csr_auto,
    two_fermion_csr,
)

from ._oracles import fixed_particle_basis, one_body_oracle, two_body_oracle


def test_fock_basis_encode_decode_roundtrip():
    """Check basis indexing used by ``ops(..., backend="scipy")``.

    Sparse Hamiltonians store rows and columns by positions in a many-electron
    basis.  This test ensures basis states can be converted to indices and back
    before matrix elements are assembled for any SciPy solver workflow.
    """
    basis = FockBasis([0b1100, 0b1010, 0b1001], norbs=4)

    assert len(basis) == 3
    for position, state in enumerate(basis.basis_int):
        assert basis.encode(state) == position
        assert basis.decode(position) == state


def test_one_body_csr_matches_independent_jordan_wigner_oracle():
    """Check one-body Hamiltonian construction in the SciPy sparse workflow.

    After ``setup_*`` defines levels, hopping, fields, and spin-orbit terms,
    ``ops(..., backend="scipy")`` lifts them into the many-electron basis.  This
    test compares that stage with an independent construction before ED/XAS/RIXS.
    """
    rng = np.random.default_rng(4)
    basis = fixed_particle_basis(4, 2)
    emat = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))

    actual = two_fermion_csr(emat, basis).toarray()
    expected = one_body_oracle(emat, basis)

    assert_allclose(actual, expected, rtol=0, atol=1e-13)


def test_one_particle_sector_reproduces_orbital_matrix():
    """Check the simplest setup-to-sparse-operator mapping.

    With one particle, the Fock-space one-body operator should be exactly the
    orbital matrix supplied by ``setup_*``.  This catches ordering mistakes
    before the resulting Hamiltonian reaches ``ed_krylov_scipy``.
    """
    basis = FockBasis([0b10, 0b01], norbs=2)
    emat = np.array([[1.2, 0.3 + 0.4j], [0.3 - 0.4j, -0.2]])

    assert_allclose(two_fermion_csr(emat, basis).toarray(), emat)


def test_four_body_dense_and_flat_sparse_paths_match_oracle():
    """Check Coulomb-operator construction in the SciPy sparse workflow.

    ``setup_*`` can provide interactions as a dense four-index tensor or as a
    flattened sparse matrix.  Both must generate the same many-electron
    interaction before initial and intermediate Hamiltonians are solved.
    """
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


def test_tolerance_drops_small_terms():
    """Check numerical pruning during sparse Hamiltonian construction.

    ``ops(..., backend="scipy")`` discards orbital matrix elements below a
    tolerance to keep operators sparse.  This test ensures only negligible
    terms are removed before the Hamiltonian is passed to a Krylov solver.
    """
    basis = fixed_particle_basis(2, 1)
    emat = np.diag([1e-11, 2e-9])

    actual = two_fermion_csr(emat, basis, tol=1e-10).toarray()

    assert actual[0, 0] == 0
    assert actual[1, 1] == pytest.approx(2e-9)


def test_basis_norbs_mismatch_raises():
    """Reject incompatible basis sectors before a transition operator is built.

    Initial and intermediate basis objects must describe the same orbital
    space when sparse one-body or transition operators are assembled.  This
    test catches an invalid setup-to-``ops`` handoff before spectroscopy starts.
    """
    left = FockBasis([0b10, 0b01], norbs=2)
    right = FockBasis([0b100], norbs=3)

    with pytest.raises(ValueError, match="same norbs"):
        two_fermion_csr(np.eye(2), left, right)


def test_sparse_u_shape_is_validated():
    """Reject a malformed sparse Coulomb tensor at the ``setup_*``/``ops`` boundary.

    Flattened sparse interactions must have one row and column for every pair
    of orbitals.  This test prevents an incorrectly shaped tensor from creating
    a wrong many-electron Hamiltonian and therefore wrong spectra.
    """
    basis = fixed_particle_basis(3, 2)

    with pytest.raises(ValueError, match="expected"):
        _four_fermion_csr_from_sparse_umat(sp.eye(8), basis)
