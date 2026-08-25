"""Numerical consistency tests for explicit/combinadic and vanilla/Numba builders."""

import importlib.util
from itertools import product

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

from edrixs.fock_basis import FockBasisSpec, build_fock_basis
from edrixs.solvers import build_op


HAS_NUMBA = importlib.util.find_spec("numba") is not None


ROUTES = [
    pytest.param(method, use_numba, id=f"{method}-{'numba' if use_numba else 'vanilla'}")
    for method, use_numba in product(("combinadic", "explicit"), (False, True))
]


def _annihilation_operators(norbs):
    """Independent full-Fock-space Jordan-Wigner annihilation matrices."""
    identity = np.eye(2, dtype=complex)
    parity = np.diag([1.0, -1.0]).astype(complex)
    annihilate = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)

    result = []
    for orbital in range(norbs):
        factors = (
            [parity] * orbital
            + [annihilate]
            + [identity] * (norbs - orbital - 1)
        )
        operator = factors[0]
        for factor in factors[1:]:
            operator = np.kron(operator, factor)
        result.append(operator)
    return result


def _oracle(emat, umat, left_spec, right_spec=None):
    """Build an independent dense Jordan-Wigner reference matrix."""
    if right_spec is None:
        right_spec = left_spec

    left = build_fock_basis(left_spec, method="explicit")
    right = build_fock_basis(right_spec, method="explicit")
    norbs = left.norbs
    annihilators = _annihilation_operators(norbs)
    full = np.zeros((2**norbs, 2**norbs), dtype=complex)

    if emat is not None:
        for iorb in range(norbs):
            for jorb in range(norbs):
                full += (
                    emat[iorb, jorb]
                    * annihilators[iorb].conj().T
                    @ annihilators[jorb]
                )

    if umat is not None:
        dense_u = umat.toarray().reshape((norbs,) * 4) if sp.issparse(umat) else umat
        for lorb, korb, jorb, iorb in zip(*np.nonzero(dense_u)):
            full += (
                dense_u[lorb, korb, jorb, iorb]
                * annihilators[lorb].conj().T
                @ annihilators[korb].conj().T
                @ annihilators[jorb]
                @ annihilators[iorb]
            )

    return full[np.ix_(left.basis_int, right.basis_int)]


def _skip_missing_numba(use_numba):
    if use_numba and not HAS_NUMBA:
        pytest.skip("numba is required for the requested JIT route")


@pytest.mark.parametrize("basis_method,use_numba", ROUTES)
def test_square_one_and_two_body_operator_matches_independent_oracle(
    basis_method, use_numba
):
    """Every SciPy construction route must produce the same physical matrix."""
    _skip_missing_numba(use_numba)
    rng = np.random.default_rng(20260817)
    raw = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    emat = (raw + raw.conj().T) / 2

    umat = np.zeros((4, 4, 4, 4), dtype=complex)
    umat[0, 1, 1, 0] = 0.7
    umat[1, 0, 0, 1] = 0.7
    umat[3, 2, 1, 0] = -0.15 + 0.11j
    umat[0, 1, 2, 3] = -0.15 - 0.11j

    spec = FockBasisSpec.from_args(4, 2)
    expected = _oracle(emat, umat, spec)
    actual = build_op(
        emat,
        umat,
        spec,
        backend="scipy",
        basis_method=basis_method,
        use_numba=use_numba,
    ).toarray()

    assert_allclose(actual, expected, rtol=0, atol=2e-13)


@pytest.mark.parametrize("basis_method,use_numba", ROUTES)
def test_rectangular_transition_operator_matches_independent_oracle(
    basis_method, use_numba
):
    """Left/right sector changes must agree for all representation/JIT routes."""
    _skip_missing_numba(use_numba)
    right = FockBasisSpec.from_args(2, 1, 2, 2)
    left = FockBasisSpec.from_args(2, 2, 2, 1)
    emat = np.zeros((4, 4), dtype=complex)
    emat[0, 2] = 1.2 - 0.3j
    emat[1, 3] = -0.4 + 0.2j

    expected = _oracle(emat, None, left, right)
    actual = build_op(
        emat,
        None,
        left,
        right,
        backend="scipy",
        basis_method=basis_method,
        use_numba=use_numba,
    ).toarray()

    assert actual.shape == expected.shape
    assert_allclose(actual, expected, rtol=0, atol=2e-13)


@pytest.mark.parametrize("basis_method,use_numba", ROUTES)
def test_dense_and_sparse_coulomb_inputs_are_numerically_identical(
    basis_method, use_numba
):
    _skip_missing_numba(use_numba)
    spec = FockBasisSpec.from_args(4, 2)
    umat = np.zeros((4, 4, 4, 4), dtype=complex)
    umat[0, 1, 1, 0] = 0.9
    umat[3, 2, 1, 0] = 0.1 - 0.05j
    sparse_u = sp.csr_matrix(umat.reshape(16, 16))

    dense_operator = build_op(
        None,
        umat,
        spec,
        backend="scipy",
        basis_method=basis_method,
        use_numba=use_numba,
    )
    sparse_operator = build_op(
        None,
        sparse_u,
        spec,
        backend="scipy",
        basis_method=basis_method,
        use_numba=use_numba,
    )

    assert_allclose(dense_operator.toarray(), sparse_operator.toarray(), atol=2e-13)
