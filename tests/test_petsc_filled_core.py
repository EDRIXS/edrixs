"""Filled-core models and centered PETSc solvers, also runnable under MPI."""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

import edrixs

PETSc = pytest.importorskip('petsc4py.PETSc')
pytest.importorskip('slepc4py.SLEPc')


def _matrix(array):
    """Distribute a small reference matrix without duplicating off-rank entries."""
    matrix = PETSc.Mat().createAIJ(array.shape, comm=PETSc.COMM_WORLD)
    matrix.setUp()
    start, end = matrix.getOwnershipRange()
    columns = np.arange(array.shape[1], dtype=PETSc.IntType)
    for row in range(start, end):
        matrix.setValues(row, columns, array[row])
    matrix.assemble()
    return matrix


def _local_values(matrix):
    start, end = matrix.getOwnershipRange()
    return matrix.getValues(
        np.arange(start, end, dtype=PETSc.IntType),
        np.arange(matrix.getSize()[1], dtype=PETSc.IntType),
    )


@pytest.mark.parametrize('method', ['explicit', 'combinadic'])
@pytest.mark.parametrize('sparse_u', [False, True])
@pytest.mark.parametrize('use_numba', [False, True])
def test_filled_core_model_build_and_shifted_ed(method, sparse_u, use_numba):
    """Core interactions survive native distributed assembly and shifted ED."""
    if use_numba:
        pytest.importorskip('numba')
    *problem, shift = edrixs.model_1v1c(
        ('p', 's'), shell_level=(0.25, -2.**20), v_noccu=1,
        slater=([0.5] * 10, [0.75] * 10), v_soc=(0.25, 0.5),
        v_othermat=np.diag(np.arange(6) / 8), sparse_U=sparse_u,
    )
    hi, hn, transitions = edrixs.get_ops(
        *problem, backend='petsc', basis_method=method, use_numba=use_numba,
        backend_kws={'assembly_chunk_cols': 2},
    )
    reference = edrixs.get_ops(*problem, backend='scipy')
    for actual, expected in zip([hi, hn, *transitions], [*reference[:2], *reference[2]]):
        start, end = actual.getOwnershipRange()
        assert_allclose(_local_values(actual), expected.toarray()[start:end], atol=1e-12)

    before = _local_values(hi)
    values, vectors = edrixs.ed(hi, num_evals=2, shift=shift)
    centered = reference[0].toarray() - shift * np.eye(hi.getSize()[0])
    expected, _ = np.linalg.eigh(centered)
    assert_allclose(values, expected[:2] + shift, rtol=0, atol=1e-9)
    residual = hi.createVecLeft()
    for value, vector in zip(values, vectors):
        hi.mult(vector, residual)
        residual.axpy(-value, vector)
        assert residual.norm() < 1e-8
    assert_array_equal(_local_values(hi), before)
    for obj in [residual, *vectors, hi, hn, *transitions]:
        obj.destroy()


@pytest.mark.parametrize('shift', [-2.**30, 2.**30])
def test_shifted_ed_accepts_zero_centered_eigenvalue(shift):
    matrix = _matrix(np.diag(np.arange(6) / 4 + shift))
    values, vectors = edrixs.ed(matrix, num_evals=2, shift=shift)
    assert_allclose(values, [shift, shift + .25], rtol=0, atol=1e-9)
    for obj in [matrix, *vectors]:
        obj.destroy()


@pytest.mark.parametrize('ntrans', [3, 5])
@pytest.mark.parametrize('skip_gs', [False, True])
@pytest.mark.parametrize('offset_i,offset_n', [
    (0., 0.), (2.**40, 2.**40), (-6 * 2.**30, -5 * 2.**30),
])
def test_centered_petsc_spectra_match_scipy(ntrans, skip_gs, offset_i, offset_n):
    """Compare complex polarization, state ordering, and absolute pole exports."""
    hi = np.diag([0., .25, .75, 1.5])
    hn = np.array([
        [2., .125j, .0625, 0], [-.125j, 2.5, .25, .125],
        [.0625, .25, 3., .125j], [0, .125, -.125j, 3.5],
    ])
    rng = np.random.default_rng(41)
    transitions = [
        (rng.integers(-8, 9, (4, 4)) + 1j * rng.integers(-8, 9, (4, 4))) / 8
        for _ in range(ntrans)
    ]
    energies = np.array([.25, 0.])
    vectors = np.eye(4, dtype=complex)[:, [1, 0]] * [1j, -1j]
    ominc = np.array([2.2, 2.7] if offset_i == offset_n else [2.25, 2.75])
    eloss = np.linspace(-.1, 1.7, 31)
    matrices = [_matrix(hi + offset_i * np.eye(4)),
                _matrix(hn + offset_n * np.eye(4))]
    ops = [_matrix(t) for t in transitions]
    states = []
    for column in vectors.T:
        vector = matrices[0].createVecRight()
        start, end = vector.getOwnershipRange()
        vector.getArray()[:] = column[start:end]
        states.append(vector)
    before = [_local_values(matrix) for matrix in matrices]
    options = dict(gamma_c=[.125, .25], temperature=3000.,
                   thin=.4, phi=.2, backend_kws={'nkryl': 4})
    xpol = [('isotropic', 0), ('linear', .3), ('left', 0), ('right', 0)]
    expected = edrixs.xas(energies, vectors, hn, transitions, ominc,
                          pol_type=xpol, backend='scipy', **options)
    actual = edrixs.xas(energies + offset_i, states, matrices[1], ops,
                        ominc + (offset_n - offset_i), pol_type=xpol, **options)
    assert_allclose(actual, expected, rtol=1e-10, atol=1e-11)

    options.update(gamma_f=np.linspace(.04, .06, len(eloss)), skip_gs=skip_gs,
                   pol_type=[('linear', .2, 'left', 0), ('right', 0, 'linear', .7)])
    expected = edrixs.rixs(
        energies, vectors, hi, hn, transitions, ominc, eloss, backend='scipy', **options,
    )
    actual, poles = edrixs.rixs(
        energies + offset_i, states, *matrices, ops, ominc + (offset_n - offset_i), eloss,
        return_poles=True, **options,
    )
    assert_allclose(actual, expected, rtol=1e-9, atol=1e-10)
    for iom, channels in enumerate(poles):
        for ip, channel in enumerate(channels):
            assert_array_equal(channel['eigval'], energies + offset_i)
            for alpha in channel['alpha']:
                assert np.all(np.abs(alpha - offset_i) < 2)
            if offset_i == 0:
                reconstructed = edrixs.get_spectra_from_poles(
                    channel, eloss, options['gamma_f'], options['temperature'])
                assert_allclose(reconstructed, actual[iom, :, ip])
    for matrix, original in zip(matrices, before):
        assert_array_equal(_local_values(matrix), original)
    for obj in [*matrices, *ops, *states]:
        obj.destroy()


def test_centered_petsc_zero_transitions():
    hi, hn = _matrix(np.diag([-100., -99.])), _matrix(np.diag([-80., -79.]))
    transitions = [_matrix(np.zeros((2, 2))) for _ in range(3)]
    state = hi.createVecRight()
    state.set(0)
    start, end = state.getOwnershipRange()
    if start == 0 and end > 0:
        state.getArray()[0] = 1
    xas = edrixs.xas([-100.], [state], hn, transitions, [20.])
    rixs, poles = edrixs.rixs(
        [-100.], [state], hi, hn, transitions, [20.], [0., 1.],
        skip_gs=True, return_poles=True,
    )
    assert_array_equal(xas, 0)
    assert_array_equal(rixs, 0)
    assert poles[0][0]['norm'] == [0]
    assert_array_equal(poles[0][0]['alpha'][0], [0])
    assert poles[0][0]['eigval'] == [-100.]
    for obj in [hi, hn, *transitions, state]:
        obj.destroy()
