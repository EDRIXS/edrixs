"""Numerical and solve-count checks for shared RIXS correction vectors.

The PETSc cases can also run with mpirun -n 2 python -m pytest on this file.
"""

import importlib
from functools import wraps

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest

import edrixs
from edrixs._solvers_helpers import _group_rixs_incoming


def test_incoming_groups_use_exact_coefficients_and_preserve_indices():
    a = np.array([1., 1j, 0.])
    near = a.copy()
    near[0] += 1e-14
    outgoing = [np.eye(3)[i % 3] for i in range(5)]
    groups = _group_rixs_incoming(list(zip([a, near, a.copy(), -a, a], outgoing)))
    assert [[i for i, _ in channels] for _, channels in groups] == [[0, 2, 4], [1], [3]]
    for _, channels in groups:
        for index, vector in channels:
            assert_array_equal(vector, outgoing[index])
    assert _group_rixs_incoming([]) == []


@pytest.fixture(params=['scipy', 'petsc', 'shifted_pminres'])
def runner(request, monkeypatch):
    """Run a small dense input problem through either backend; count real solves."""
    if request.param == 'scipy':
        backend = importlib.import_module('edrixs.scipy_backend.scipy_backend')
        original = backend.gmres
        counts = []

        @wraps(original)
        def counted(*args, **kwargs):
            counts.append(1)
            return original(*args, **kwargs)

        monkeypatch.setattr(backend, 'gmres', counted)
        return backend.rixs_scipy, lambda: len(counts)

    PETSc = pytest.importorskip('petsc4py.PETSc')
    options = PETSc.Options()
    saved = options.getAll()
    options['pc_type'] = 'none'
    options['ksp_rtol'] = 0.0

    def restore_options():
        for key in ('pc_type', 'ksp_rtol'):
            if key in saved:
                options[key] = saved[key]
            else:
                del options[key]

    request.addfinalizer(restore_options)
    PETSc.Log.begin()
    stage = PETSc.Log.Stage('test_rixs_reuse')
    event = PETSc.Log.Event('KSPSolve')
    initial_count = event.getPerfInfo(stage.id)['count']

    def run(energies, vectors, hi, hn, transitions, ominc, eloss, **kwargs):
        if request.param == 'petsc':
            kwargs['backend_kws'] = dict(kwargs.get('backend_kws', {}),
                                         ksp_type='gmres')
        objects = []

        def matrix(array):
            mat = PETSc.Mat().createAIJ(array.shape, comm=PETSc.COMM_WORLD)
            objects.append(mat)
            mat.setUp()
            start, end = mat.getOwnershipRange()
            columns = np.arange(array.shape[1], dtype=PETSc.IntType)
            for row in range(start, end):
                mat.setValues(row, columns, array[row])
            mat.assemble()
            return mat

        try:
            initial, intermediate = matrix(hi), matrix(hn)
            ops = [matrix(t) for t in transitions]
            states = []
            for column in vectors.T:
                state = initial.createVecRight()
                objects.append(state)
                start, end = state.getOwnershipRange()
                state.getArray()[:] = column[start:end]
                states.append(state)
            originals = [state.getArray(readonly=True).copy() for state in states]
            stage.push()
            try:
                result = edrixs.rixs(energies, states, initial, intermediate, ops,
                                     ominc, eloss, backend='petsc', **kwargs)
            finally:
                stage.pop()
            for state, original in zip(states, originals):
                assert_array_equal(state.getArray(readonly=True), original)
            return result
        finally:
            for obj in reversed(objects):
                obj.destroy()

    return run, lambda: event.getPerfInfo(stage.id)['count'] - initial_count


def problem(ntrans, degenerate=False):
    energies = np.array([0. if degenerate else .25, 0.])
    hi = np.diag([0., energies[0], .75, 1.5])
    hn = np.array([[2., .125j, .0625, 0], [-.125j, 2.5, .25, .125],
                   [.0625, .25, 3., .125j], [0, .125, -.125j, 3.5]])
    rng = np.random.default_rng(42)
    transitions = [(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
                   for _ in range(ntrans)]
    vectors = np.eye(4, dtype=complex)[:, [1, 0]] * [1j, -1j]
    # Same incident energy with different broadening must not share a solve.
    return energies, vectors, hi, hn, transitions, [2.2, 2.2], np.linspace(-.1, 1.7, 21)


@pytest.mark.parametrize('ntrans', [3, 5])
@pytest.mark.parametrize('skip_gs', [False, True])
@pytest.mark.parametrize('degenerate', [False, True])
def test_reuse_matches_single_channels_and_dense_reference(runner, ntrans, skip_gs, degenerate):
    run, solve_count = runner
    args = problem(ntrans, degenerate)
    channels = [('left', 0., 'right', 0.), ('linear', .2, 'linear', .7),
                ('left', 0., 'linear', .4), ('left', 0., 'right', 0.)]
    options = dict(gamma_c=[.125, .25], gamma_f=.05, temperature=3000.,
                   thin=.4, thout=.6, phi=.2, skip_gs=skip_gs)
    solver_options = {'nkryl': 4, 'linsys_tol': 1e-12}
    actual, poles = run(*args, pol_type=channels, return_poles=True,
                        backend_kws=solver_options, **options)
    # Two energies x two initial states x two distinct incoming polarizations.
    assert solve_count() == 8
    dense = edrixs.rixs(*args, pol_type=channels, backend='dense', **options)
    assert_allclose(actual, dense, rtol=1e-8, atol=1e-9)
    for ip, channel in enumerate(channels):
        separate = run(*args, pol_type=[channel], backend_kws=solver_options, **options)
        assert_allclose(actual[:, :, ip], separate[:, :, 0], rtol=1e-10, atol=1e-11)
        for iom in range(2):
            record = poles[iom][ip]
            assert_array_equal(record['eigval'], args[0])
            rebuilt = edrixs.get_spectra_from_poles(record, args[-1], .05, 3000.)
            assert_allclose(actual[iom, :, ip], rebuilt, rtol=1e-11, atol=1e-11)
    assert solve_count() == 8 + 16
    # Duplicate outputs must not alias mutable pole arrays.
    assert poles[0][0]['alpha'][0] is not poles[0][3]['alpha'][0]


def test_zero_incoming_vectors_bypass_solves_for_all_channels(runner):
    run, solve_count = runner
    args = list(problem(3))
    args[4] = [np.zeros((4, 4)) for _ in range(3)]
    actual, poles = run(*args, pol_type=[('left', 0, 'left', 0), ('left', 0, 'right', 0)],
                        return_poles=True, skip_gs=True)
    assert solve_count() == 0
    assert_array_equal(actual, 0)
    for channels in poles:
        for record in channels:
            assert record['norm'] == [0., 0.]
            assert record['npoles'] == [1, 1]
            for alpha, beta in zip(record['alpha'], record['beta']):
                assert_array_equal(alpha, [0.])
                assert beta.size == 0


def test_failed_solve_does_not_continue_to_outgoing_channels(runner):
    run, _ = runner
    with pytest.raises(RuntimeError, match='did not converge'):
        run(*problem(3), pol_type=[('left', 0, 'left', 0), ('left', 0, 'right', 0)],
            backend_kws={'linsys_maxiter': 1, 'linsys_tol': 1e-30, 'nkryl': 4})
