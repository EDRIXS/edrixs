"""Unit and numerical consistency tests for the PETSc backend."""

import importlib.util

import numpy as np
import pytest
from numpy.testing import assert_allclose

import edrixs.petsc_backend.petsc_backend as backend
from edrixs.fock_basis import FockBasisSpec
from edrixs.petsc_backend.lanczos import lanczos_tridiagonal as lanczos_petsc
from edrixs.scipy_backend.krylov import lanczos_tridiagonal as lanczos_scipy
from edrixs.solvers import build_op, ed


HAS_PETSC = importlib.util.find_spec("petsc4py") is not None
HAS_SLEPC = importlib.util.find_spec("slepc4py") is not None
requires_petsc = pytest.mark.skipif(
    not HAS_PETSC,
    reason="petsc4py is required for PETSc numeric tests",
)
requires_petsc_slepc = pytest.mark.skipif(
    not (HAS_PETSC and HAS_SLEPC),
    reason="petsc4py and slepc4py are required for PETSc eigensolver tests",
)


# -----------------------------------------------------------------------------
# Contract / logic tests (no PETSc runtime required)
# -----------------------------------------------------------------------------


def test_owns_operator_returns_false_for_non_petsc_objects():
    """Recognition must never raise when given ordinary Python/NumPy objects."""
    assert backend.owns_operator_petsc(np.eye(2)) is False
    assert backend.owns_operator_petsc(object()) is False


def test_backend_kws_validation():
    assert backend._backend_kws(None) == {}
    assert backend._backend_kws({"nkryl": 50}) == {"nkryl": 50}
    with pytest.raises(TypeError, match="mapping"):
        backend._backend_kws([("nkryl", 50)])


def test_backend_kws_returns_independent_copy():
    original = {"nkryl": 10}
    copied = backend._backend_kws(original)
    copied["nkryl"] = 999
    assert original["nkryl"] == 10


@pytest.mark.parametrize("num_evals", [0, -3])
def test_ed_petsc_rejects_nonpositive_num_evals(num_evals):
    if HAS_PETSC and HAS_SLEPC:
        pytest.skip("covered by numeric tests when PETSc/SLEPc are present")
    with pytest.raises((ValueError, ImportError)):
        backend.ed_petsc(object(), num_evals=num_evals)


def test_rixs_petsc_rejects_unknown_options_without_petsc():
    if HAS_PETSC:
        pytest.skip("import guard is only meaningful without petsc4py")
    with pytest.raises(ImportError):
        backend.rixs_petsc(
            np.array([0.0]),
            [object()],
            object(),
            object(),
            [object()] * 3,
            np.array([1.0]),
            np.array([0.0]),
            backend_kws={"not_a_real_option": 1},
        )


def test_public_symbols_are_exported():
    for name in (
        "owns_operator_petsc",
        "build_op_petsc",
        "ed_petsc",
        "xas_petsc",
        "rixs_petsc",
    ):
        assert name in backend.__all__
        assert hasattr(backend, name)


# -----------------------------------------------------------------------------
# Numeric PETSc tests
# -----------------------------------------------------------------------------


def _petsc_dense_matrix(array):
    from petsc4py import PETSc

    array = np.asarray(array, dtype=complex)
    size = array.shape[0]
    matrix = PETSc.Mat().createAIJ((size, size), comm=PETSc.COMM_SELF)
    matrix.setUp()
    columns = np.arange(size, dtype=PETSc.IntType)
    for row in range(size):
        matrix.setValues(row, columns, array[row])
    matrix.assemblyBegin()
    matrix.assemblyEnd()
    return matrix


def _petsc_vector(matrix, values):
    from petsc4py import PETSc

    vector = matrix.createVecRight()
    indices = np.arange(len(values), dtype=PETSc.IntType)
    vector.setValues(indices, np.asarray(values, dtype=complex))
    vector.assemblyBegin()
    vector.assemblyEnd()
    return vector


@requires_petsc
def test_petsc_lanczos_matches_scipy_lanczos_and_preserves_seed():
    """Backend-specific Lanczos implementations must produce the same projection."""
    rng = np.random.default_rng(91)
    raw = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6))
    hermitian = (raw + raw.conj().T) / 2
    seed = rng.normal(size=6) + 1j * rng.normal(size=6)

    matrix = _petsc_dense_matrix(hermitian)
    vector = _petsc_vector(matrix, seed)
    vector_before = vector.copy()

    alpha_p, beta_p, norm_p = lanczos_petsc(matrix, vector, nkryl=20)
    alpha_s, beta_s, norm_s = lanczos_scipy(hermitian, seed, m=20)

    assert_allclose(alpha_p, alpha_s, rtol=0, atol=2e-12)
    assert_allclose(beta_p, beta_s, rtol=0, atol=2e-12)
    assert norm_p == pytest.approx(norm_s, abs=2e-12)
    vector.axpy(-1.0, vector_before)
    assert vector.norm() == pytest.approx(0.0, abs=1e-14)


@requires_petsc_slepc
def test_public_ed_scipy_and_petsc_agree_on_constructed_operator():
    """The two solver backends agree after constructing the same Fock operator."""
    from petsc4py import PETSc

    rng = np.random.default_rng(123)
    raw = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    emat = (raw + raw.conj().T) / 2
    umat = np.zeros((4, 4, 4, 4), dtype=complex)
    umat[0, 1, 1, 0] = 0.35
    umat[1, 0, 0, 1] = 0.35
    spec = FockBasisSpec.from_args(4, 2)

    scipy_matrix = build_op(
        emat,
        umat,
        spec,
        backend="scipy",
        basis_method="combinadic",
        use_numba=False,
    )
    petsc_matrix = build_op(
        emat,
        umat,
        spec,
        backend="petsc",
        basis_method="combinadic",
        use_numba=False,
        backend_kws={"comm": PETSc.COMM_SELF},
    )

    eval_s, _ = ed(scipy_matrix, num_evals=2, backend="scipy")
    eval_p, _ = ed(petsc_matrix, num_evals=2, backend="petsc")

    assert_allclose(eval_p, eval_s, rtol=0, atol=1e-8)


@pytest.fixture
def complex_hermitian_petsc_mat():
    if not HAS_PETSC:
        pytest.skip("petsc4py is required")
    rng = np.random.default_rng(12345)
    dense = rng.standard_normal((6, 6))
    dense = dense + 1j * rng.standard_normal((6, 6))
    hermitian = (dense + dense.conj().T) / 2
    return _petsc_dense_matrix(hermitian), np.linalg.eigvalsh(hermitian)


@requires_petsc_slepc
def test_ed_petsc_returns_lowest_sorted_eigenpairs(complex_hermitian_petsc_mat):
    mat, expected_eigenvalues = complex_hermitian_petsc_mat
    evals, evecs = backend.ed_petsc(mat, num_evals=3)

    assert len(evals) == 3
    assert len(evecs) == 3
    assert_allclose(evals, expected_eigenvalues[:3], atol=1e-8)

    residual = mat.getVecLeft()
    for value, vector in zip(evals, evecs):
        mat.mult(vector, residual)
        residual.axpy(-value, vector)
        assert residual.norm() < 1e-6


@requires_petsc_slepc
def test_ed_petsc_rejects_too_many_requested(complex_hermitian_petsc_mat):
    mat, expected_eigenvalues = complex_hermitian_petsc_mat
    with pytest.raises((RuntimeError, ValueError)):
        backend.ed_petsc(mat, num_evals=len(expected_eigenvalues) + 5)
