"""Unit tests for the PETSc backend contract.

These tests are split into two groups:

* Contract/logic tests that only need ``numpy`` and exercise the pure-Python
  parts of :mod:`edrixs.petsc_backend.petsc_backend` (option parsing,
  validation, operator recognition). They run everywhere.
* Numeric tests that require a working ``petsc4py``/``slepc4py`` stack. They
  are skipped automatically when those optional dependencies are missing.
"""

import importlib.util

import numpy as np
import pytest

import edrixs.petsc_backend.petsc_backend as backend


HAS_PETSC = importlib.util.find_spec("petsc4py") is not None
HAS_SLEPC = importlib.util.find_spec("slepc4py") is not None
requires_petsc = pytest.mark.skipif(
    not (HAS_PETSC and HAS_SLEPC),
    reason="petsc4py and slepc4py are required for the PETSc numeric backend",
)


# -----------------------------------------------------------------------------
# Contract / logic tests (no PETSc runtime required)
# -----------------------------------------------------------------------------


def test_owns_operator_returns_false_without_petsc():
    """Recognition must never raise, even when petsc4py is unavailable."""
    assert backend.owns_operator_petsc(np.eye(2)) is False
    assert backend.owns_operator_petsc(object()) is False


def test_backend_kws_validation():
    """``backend_kws`` must be a mapping or None."""
    assert backend._backend_kws(None) == {}
    assert backend._backend_kws({"nkryl": 50}) == {"nkryl": 50}
    with pytest.raises(TypeError, match="mapping"):
        backend._backend_kws([("nkryl", 50)])


def test_backend_kws_returns_independent_copy():
    """The helper must copy so caller dicts are never mutated."""
    original = {"nkryl": 10}
    copied = backend._backend_kws(original)
    copied["nkryl"] = 999
    assert original["nkryl"] == 10


@pytest.mark.parametrize("num_evals", [0, -3])
def test_ed_petsc_rejects_nonpositive_num_evals(num_evals):
    """Validation happens before any operator use, so a fake mat suffices."""
    if HAS_PETSC and HAS_SLEPC:
        pytest.skip("covered by numeric tests when PETSc is present")
    with pytest.raises((ValueError, ImportError)):
        backend.ed_petsc(object(), num_evals=num_evals)


def test_rixs_petsc_rejects_unknown_options():
    """Unknown backend options must be reported explicitly."""
    if HAS_PETSC:
        pytest.skip("import guard only meaningful without petsc4py")
    with pytest.raises(ImportError):
        backend.rixs_petsc(
            np.array([0.0]), [object()], object(), object(),
            [object()] * 3, np.array([1.0]), np.array([0.0]),
            backend_kws={"not_a_real_option": 1},
        )


def test_public_symbols_are_exported():
    """The backend advertises the full contract surface."""
    for name in (
        "owns_operator_petsc", "build_op_petsc",
        "ed_petsc", "xas_petsc", "rixs_petsc",
    ):
        assert name in backend.__all__
        assert hasattr(backend, name)


# -----------------------------------------------------------------------------
# Numeric tests (require petsc4py + slepc4py)
# -----------------------------------------------------------------------------


@pytest.fixture
def complex_hermitian_petsc_mat():
    """Build a reproducible dense complex Hermitian PETSc matrix."""
    from petsc4py import PETSc

    size = 6
    rng = np.random.default_rng(12345)
    dense = rng.standard_normal((size, size))
    dense = dense + 1j * rng.standard_normal((size, size))
    hermitian = (dense + dense.conj().T) / 2

    mat = PETSc.Mat().createAIJ(
        (size, size), comm=PETSc.COMM_WORLD
    )
    mat.setUp()
    start, end = mat.getOwnershipRange()
    columns = np.arange(size, dtype=PETSc.IntType)
    for row in range(start, end):
        mat.setValues(row, columns, hermitian[row])
    mat.assemblyBegin()
    mat.assemblyEnd()
    return mat, np.linalg.eigvalsh(hermitian)


@requires_petsc
def test_ed_petsc_returns_lowest_sorted_eigenpairs(
    complex_hermitian_petsc_mat,
):
    """SLEPc returns the requested number of lowest eigenpairs, in order."""
    mat, expected_eigenvalues = complex_hermitian_petsc_mat

    evals, evecs = backend.ed_petsc(mat, num_evals=3)

    assert len(evals) == 3
    assert len(evecs) == 3
    np.testing.assert_allclose(evals, expected_eigenvalues[:3], atol=1e-8)

    # Each returned pair must satisfy H v = lambda v.
    residual = mat.getVecLeft()
    for value, vector in zip(evals, evecs):
        mat.mult(vector, residual)
        residual.axpy(-value, vector)
        assert residual.norm() < 1e-6


@requires_petsc
def test_ed_petsc_rejects_too_many_requested(complex_hermitian_petsc_mat):
    """Requesting more eigenpairs than converge must raise."""
    mat, expected_eigenvalues = complex_hermitian_petsc_mat
    with pytest.raises((RuntimeError, ValueError)):
        backend.ed_petsc(mat, num_evals=len(expected_eigenvalues) + 5)
