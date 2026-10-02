"""Explicit filled-core models and their native Fortran disk reduction."""

from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose, assert_array_equal

import edrixs
from edrixs.fortran_backend import fortran_backend as fort
from edrixs.fortran_backend.iostream_fortran import write_emat, write_umat


def _integrals():
    rng = np.random.default_rng(214)
    h = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6))
    h += h.conj().T
    u = rng.normal(size=(6,) * 4) + 1j * rng.normal(size=(6,) * 4)
    u += u.transpose(3, 2, 1, 0).conj()
    return h, u


@pytest.mark.parametrize('nocc', [0, 1, 2, 3, 4])
@pytest.mark.parametrize('sparse', [False, True])
def test_reduction_matches_explicit_many_body_hamiltonian(nocc, sparse):
    """Exercise all contractions with complex, non-antisymmetrized integrals."""
    h, u = _integrals()
    if sparse:
        u = sp.csr_matrix(u.reshape(36, 36))
    h_before, u_before = h.copy(), u.copy()
    effective, valence_u, constant = fort._reduce_initial_integrals(h, u, 4)
    full = edrixs.build_op(
        h, u, edrixs.FockBasisSpec.from_args(4, nocc, 2, 2), backend='dense')
    reduced = edrixs.build_op(
        effective, valence_u, edrixs.FockBasisSpec.from_args(4, nocc), backend='dense')
    assert_allclose(full, reduced + constant * np.eye(len(reduced)), atol=2e-13)
    assert abs(constant) > 0.1
    assert np.max(np.abs(effective - h[:4, :4])) > 0.1
    assert_array_equal(h, h_before)
    if sparse:
        assert (u != u_before).nnz == 0
        assert sp.issparse(valence_u)
    else:
        assert_array_equal(u, u_before)


def test_sparse_staging_without_densification_and_without_mutation(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    h, u = _integrals()
    h_before, u_before = h.copy(), u.copy()
    effective, valence_u, constant = fort._reduce_initial_integrals(h, u, 4)
    write_emat(effective + constant / 2 * np.eye(4), 'expected_h.in')
    write_umat(valence_u, 'expected_u.in')
    sparse = sp.csr_matrix(u.reshape(36, 36))
    sparse_before = sparse.copy()

    def forbidden(*args, **kwargs):
        raise AssertionError('sparse Coulomb tensor was densified')

    monkeypatch.setattr(sp.csr_matrix, 'toarray', forbidden)
    monkeypatch.setattr(sp.coo_matrix, 'toarray', forbidden)
    fort._write_initial_integrals(h, sparse, 4, 2)
    # Floating-point summation order can differ between dense and sparse.
    for actual, expected in [('hopping_i.in', 'expected_h.in'),
                             ('coulomb_i.in', 'expected_u.in')]:
        assert_allclose(np.loadtxt(actual, skiprows=1),
                        np.loadtxt(expected, skiprows=1), atol=1e-13)
    assert_array_equal(h, h_before)
    assert_array_equal(u, u_before)
    assert (sparse != sparse_before).nnz == 0


@pytest.mark.parametrize('method', ['combinadic', 'explicit'])
def test_native_basis_encoding_preserves_restricted_sectors(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)
    initial = edrixs.FockBasisSpec.from_args(2, 1, 2, 0, 2, 2)
    intermediate = edrixs.FockBasisSpec.from_args(2, 1, 2, 1, 2, 1)
    h, u = _integrals()
    fort.write_problem(
        h, u, edrixs.build_fock_basis(initial, method=method),
        h, u, edrixs.build_fock_basis(intermediate, method=method),
        np.zeros((3, 6, 6)),
    )
    assert Path('fock_i.in').read_text().split() == ['2', '1', '2']
    assert Path('fock_n.in').read_text().split() == ['4', '5', '6', '9', '10']
    assert Path('fock_f.in').read_bytes() == Path('fock_i.in').read_bytes()
    # The intermediate files keep full orbital dimensions and coefficients.
    write_emat(h, 'expected_h.in')
    write_umat(u, 'expected_u.in')
    assert Path('hopping_n.in').read_bytes() == Path('expected_h.in').read_bytes()
    assert Path('coulomb_n.in').read_bytes() == Path('expected_u.in').read_bytes()


def test_zero_valence_native_limit(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    problem = edrixs.model_1v1c(('s', 's'), v_noccu=0, shell_level=(0, -2))
    full = edrixs.get_ops(*problem[:7], backend='dense')[0]
    assert_allclose(full, [[-4]])
    with pytest.raises(ValueError, match='nonzero core energy.*zero-valence'):
        fort.write_problem(*problem[:7])
    # No constant is needed for this sector, so it is representable.
    fort.write_problem(*edrixs.model_1v1c(('s', 's'), v_noccu=0)[:7])
    assert Path('fock_i.in').read_text().split() == ['1', '0']


@pytest.mark.parametrize('model,kwargs,nv,nc,nocc', [
    (edrixs.model_1v1c, dict(shell_name=('s', 'p'), v_noccu=1,
                             shell_level=(0.3, -2)), 2, 6, 1),
    (edrixs.model_2v1c, dict(shell_name=('s', 's', 'p'), v_tot_noccu=1,
                             shell_level=(0.3, 0.4, -2)), 4, 6, 1),
    (edrixs.model_siam, dict(shell_name=('s', 'p'), nbath=1, v_noccu=1,
                             c_level=-2), 4, 6, 1),
    (edrixs.model_siam_2d1p, dict(impurity_levels=np.zeros(5), bath_levels=np.zeros(5),
                                  hyb=np.zeros(5), Delta=2, nd=9), 20, 6, 19),
])
def test_model_full_shapes_core_blocks_and_sparse_agreement(model, kwargs, nv, nc, nocc):
    # Nonzero entries throughout the Slater list retain core-valence/core-core U.
    slater = ([0.3] * 30, [0.4] * 30)
    dense = model(**kwargs, slater=slater, c_soc=0.2)
    sparse = model(**kwargs, slater=slater, c_soc=0.2, sparse_U=True)
    n = nv + nc
    assert dense[0].shape == dense[3].shape == (n, n)
    assert dense[2].shapes == ((nv, nocc), (nc, nc))
    assert dense[5].shapes == ((nv, nocc + 1), (nc, nc - 1))
    assert_allclose(dense[0][nv:, nv:], dense[3][nv:, nv:])
    assert np.max(np.abs(dense[0][nv:, nv:])) > 0
    assert np.max(np.abs(dense[1][:nv, nv:, nv:, :nv])) > 0
    assert np.max(np.abs(dense[1][nv:, nv:, nv:, nv:])) > 0
    for i in (1, 4):
        assert dense[i].shape == (n,) * 4
        assert sparse[i].shape == (n*n, n*n)
        assert_allclose(sparse[i].toarray(), dense[i].reshape(n*n, n*n), atol=1e-10)
