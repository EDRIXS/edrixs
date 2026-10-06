"""Core-level model offsets and non-mutating ED centering."""

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose, assert_array_equal
from scipy.sparse.linalg import aslinearoperator

import edrixs
from edrixs import solvers


@pytest.mark.parametrize('model,kwargs,nc,core_level', [
    (edrixs.model_1v1c, dict(shell_name=('s', 'p'), shell_level=(0.2, 100)), 6, 100),
    (edrixs.model_2v1c, dict(shell_name=('s', 's', 'p'),
                             shell_level=(0.2, 0.3, -100)), 6, -100),
    (edrixs.model_siam, dict(shell_name=('s', 'p'), nbath=1, c_level=75), 6, 75),
])
@pytest.mark.parametrize('sparse', [False, True])
def test_models_append_shell_level_shift(model, kwargs, nc, core_level, sparse):
    output = model(**kwargs, slater=([0.4] * 30, [0.5] * 30),
                   c_soc=0.2, sparse_U=sparse)
    assert len(output) == 8
    *problem, shift = output
    assert isinstance(shift, float)
    assert shift == nc * core_level
    # SOC and nonzero core interactions do not change this simple estimate.
    assert_allclose(np.trace(problem[0][-nc:, -nc:]), shift)
    assert_allclose(problem[0][-nc:, -nc:], problem[3][-nc:, -nc:])
    if sp.issparse(problem[1]):
        assert problem[1].nnz > 0
    else:
        assert np.any(problem[1][-nc:, -nc:, -nc:, -nc:])


@pytest.mark.parametrize('model,args', [
    (edrixs.model_1v1c, (('s', 's'),)),
    (edrixs.model_2v1c, (('s', 's', 'p'),)),
    (edrixs.model_siam, (('s', 'p'), 1)),
])
def test_unspecified_core_level_has_zero_shift(model, args):
    assert model(*args)[-1] == 0.0


def _operator(matrix, kind):
    if kind == 'dense':
        return matrix.copy()
    if kind == 'sparse':
        return sp.csr_matrix(matrix)
    return aslinearoperator(matrix)


@pytest.mark.parametrize('kind,backend', [
    ('dense', 'dense'), ('dense', 'scipy'), ('sparse', 'scipy'), ('linear', 'scipy'),
])
@pytest.mark.parametrize('shift', [-1e6, 1e6])
@pytest.mark.parametrize('direct', [False, True])
def test_backend_centers_before_diagonalization_and_restores_energies(
        kind, backend, shift, direct, monkeypatch):
    base = np.array([[1, 0.3j], [-0.3j, 2]], dtype=complex)
    full = base + shift * np.eye(2)
    operator = _operator(full, kind)
    module = getattr(solvers, backend + '_backend')
    vectors = np.eye(2)

    def solver(centered, initial=None, **kwargs):
        if backend == 'dense':
            assert kwargs == {'subset_by_index': None}
        else:
            assert initial.shape == (2, 2)
            assert kwargs['largest'] is False
        assert centered is not operator
        assert_allclose(centered @ vectors, base)
        if kind == 'linear':
            x = np.array([0.3j, 0.7])
            block = np.column_stack([x, x.conj()])
            assert_allclose(centered.matvec(x), base @ x)
            assert_allclose(centered.rmatvec(x), base.conj().T @ x)
            assert_allclose(centered.matmat(block), base @ block)
            assert_allclose(centered.rmatmat(block), base.conj().T @ block)
        return np.array([1., 2.]), vectors

    if backend == 'dense':
        monkeypatch.setattr(module.scipy.linalg, 'eigh', solver)
    else:
        monkeypatch.setattr(module, 'lobpcg', solver)
    if kind == 'sparse':
        convert = module.aslinearoperator

        def check_sparse(matrix):
            assert sp.issparse(matrix)
            return convert(matrix)

        monkeypatch.setattr(module, 'aslinearoperator', check_sparse)

        def forbidden(*args, **kwargs):
            raise AssertionError('sparse Hamiltonian was densified')
        monkeypatch.setattr(sp.csr_matrix, 'toarray', forbidden)
    if direct:
        values, actual_vectors = getattr(module, 'ed_' + backend)(
            operator, num_evals=2, shift=shift)
    else:
        values, actual_vectors = edrixs.ed(
            operator, num_evals=2, shift=shift, backend=backend)
    assert_allclose(values, [shift + 1, shift + 2])
    assert_array_equal(actual_vectors, vectors)
    assert_array_equal(operator @ np.eye(2), full)


@pytest.mark.parametrize('kind', ['dense', 'sparse', 'linear'])
@pytest.mark.parametrize('shift', [-1e6, 1e6])
def test_shifted_ed_eigenpairs_match_original_hamiltonian(kind, shift):
    rng = np.random.default_rng(173)
    raw = rng.normal(size=(20, 20)) + 1j * rng.normal(size=(20, 20))
    base = (raw + raw.conj().T) / 2
    full = base + shift * np.eye(20)
    operator = _operator(full, kind)
    kws = {} if kind == 'dense' else {'seed': 17, 'maxiter': 200}
    energies, vectors = edrixs.ed(operator, num_evals=2, shift=shift, backend_kws=kws)
    expected_e, expected_v = np.linalg.eigh(full)
    assert_allclose(energies, expected_e[:2], atol=2e-9, rtol=0)
    assert_allclose(full @ vectors, vectors * energies, atol=2e-9, rtol=0)
    assert_allclose(vectors @ vectors.conj().T,
                    expected_v[:, :2] @ expected_v[:, :2].conj().T, atol=2e-8)
    assert_array_equal(operator @ np.eye(20), full)


@pytest.mark.parametrize('backend', ['dense', 'scipy', 'fortran', 'petsc'])
@pytest.mark.parametrize('shift', [0, -1e6, 1e6])
def test_public_ed_passes_shift_and_hamiltonian_to_backend(backend, shift, monkeypatch):
    operator = object()
    result = (np.array([2.]), object())

    def solver(actual, **kwargs):
        assert actual is operator
        assert kwargs == {'num_evals': 1, 'shift': shift, 'backend_kws': None}
        return result

    monkeypatch.setattr(getattr(solvers, backend + '_backend'), 'ed_' + backend, solver)
    assert edrixs.ed(operator, backend=backend, shift=shift) is result


@pytest.mark.parametrize('backend', ['petsc'])
@pytest.mark.parametrize('direct', [False, True])
def test_nonzero_shift_rejected_for_unsupported_backend(backend, direct):
    with pytest.raises(ValueError, match='only by dense and scipy'):
        if direct:
            module = getattr(solvers, backend + '_backend')
            getattr(module, 'ed_' + backend)(object(), shift=2)
        else:
            edrixs.ed(object(), shift=2, backend=backend)


@pytest.mark.parametrize('backend', ['dense', 'scipy'])
def test_spectra_keep_original_reference_after_shifted_ed(backend):
    *problem, shift = edrixs.model_1v1c(
        ('s', 'p'), shell_level=(0.1, 200), c_soc=0.2,
        v_othermat=np.diag([-0.2, 0.2]), slater=([0.3] * 6, [0.4] * 6))
    hi, hn, trans = edrixs.get_ops(*problem, backend=backend)
    reference_e, reference_v = edrixs.ed(hi, num_evals=2, backend=backend)
    actual_e, actual_v = edrixs.ed(hi, num_evals=2, shift=shift, backend=backend)
    assert_allclose(actual_e, reference_e, atol=1e-10)
    for energies, vectors in [(reference_e, reference_v), (actual_e, actual_v)]:
        xas = edrixs.xas(energies, vectors, hn, trans, [-201, -200, -199],
                         backend=backend, pol_type=[('isotropic', 0)])
        rixs = edrixs.rixs(energies, vectors, hi, hn, trans, [-200], [0, 0.2, 0.4],
                           backend=backend)
        if energies is reference_e:
            expected_x, expected_r = xas, rixs
        else:
            assert_allclose(xas, expected_x, atol=1e-9)
            assert_allclose(rixs, expected_r, atol=1e-9)


@pytest.mark.parametrize('direct', [False, True])
@pytest.mark.parametrize('shift', [0, -2, 2])
def test_fortran_ignores_shift_and_continues(direct, shift, tmp_path, monkeypatch):
    from contextlib import nullcontext
    import warnings
    from edrixs.fortran_backend import fortran_backend as backend

    monkeypatch.chdir(tmp_path)
    calls = []

    def fake_solver(name):
        assert name == 'ed_fsolver'

        def run(fcomm, rank, size):
            calls.append(name)
            (tmp_path / 'eigvals.dat').write_text('1 -0.25\n')

        return run

    monkeypatch.setattr(backend, '_solver', fake_solver)
    problem = edrixs.model_1v1c(('s', 's'), v_noccu=1)
    hi, _, _ = edrixs.get_ops(*problem[:7], backend='fortran')
    expected_warning = (
        pytest.warns(UserWarning, match='Fortran ED ignores shift')
        if shift != 0 else nullcontext()
    )
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        with expected_warning:
            if direct:
                energies, vectors = backend.ed_fortran(hi, shift=shift)
            else:
                energies, vectors = edrixs.ed(hi, shift=shift, backend='fortran')

    assert calls == ['ed_fsolver']
    assert_array_equal(energies, [-0.25])
    assert isinstance(vectors, backend.FortranDiskOperator)
