"""Integration tests from model metadata through all SciPy construction routes."""

import importlib.util
from itertools import product

import numpy as np
import pytest
from numpy.testing import assert_allclose

from edrixs.fock_basis import FockBasisSpec
from edrixs.models import model_1v1c, model_2v1c, model_siam
from edrixs.solvers import get_ops


HAS_NUMBA = importlib.util.find_spec("numba") is not None


# Public SciPy construction matrix: basis representation x JIT choice.
ALL_ROUTES = [
    pytest.param(
        basis_method,
        use_numba,
        id=f"scipy-{basis_method}-{'numba' if use_numba else 'vanilla'}",
    )
    for basis_method, use_numba in product(
        ("combinadic", "explicit"), (False, True)
    )
]


def _small_1v1c_problem():
    """Return a small nontrivial model shared by all route comparisons."""
    return model_1v1c(
        ("p", "s"),
        shell_level=(0.2, -4.0),
        v_soc=(0.05, 0.08),
        v_noccu=1,
        slater=([0.7], [0.9]),
        v_cfmat=np.diag(np.linspace(-0.12, 0.12, 6)),
        sparse_U=False,
    )


def _skip_unavailable(use_numba):
    """Skip only an explicitly requested optional Numba route when unavailable."""
    if use_numba and not HAS_NUMBA:
        pytest.skip("numba is required for the requested JIT route")


def test_all_model_functions_return_compact_basis_metadata():
    """No model function should materialize the many-body basis anymore."""
    problems = [
        model_1v1c(("p", "s"), v_noccu=1),
        model_2v1c(("s", "s", "p"), v_tot_noccu=1, trans_to_which=1),
        model_siam(("s", "p"), 1, v_noccu=1),
    ]

    expected_shapes = [
        (((6, 1), (2, 2)), ((6, 2), (2, 1))),
        (((4, 1), (6, 6)), ((4, 2), (6, 5))),
        (((4, 1), (6, 6)), ((4, 2), (6, 5))),
    ]

    for problem, (initial_shapes, intermediate_shapes) in zip(problems, expected_shapes):
        basis_i, basis_n = problem[2], problem[5]
        assert isinstance(basis_i, FockBasisSpec)
        assert isinstance(basis_n, FockBasisSpec)
        assert basis_i.shapes == initial_shapes
        assert basis_n.shapes == intermediate_shapes
        assert not hasattr(basis_i, "basis_int")
        assert not hasattr(basis_n, "basis_int")


@pytest.mark.integration
@pytest.mark.parametrize("basis_method,use_numba", ALL_ROUTES)
def test_model_to_get_ops_all_scipy_routes_are_numerically_consistent(
    basis_method, use_numba
):
    """Compare all four public SciPy operator-construction routes numerically.

    The explicit, non-Numba SciPy route is the common reference. Every
    combinadic/explicit x vanilla/Numba SciPy route must produce the same
    initial Hamiltonian, intermediate Hamiltonian, and transition operators
    from the same compact model metadata.
    """
    _skip_unavailable(use_numba)
    problem = _small_1v1c_problem()

    reference_i, reference_n, reference_t = get_ops(
        *problem,
        backend="scipy",
        basis_method="explicit",
        use_numba=False,
    )
    actual_i, actual_n, actual_t = get_ops(
        *problem,
        backend="scipy",
        basis_method=basis_method,
        use_numba=use_numba,
    )

    assert_allclose(actual_i.toarray(), reference_i.toarray(), rtol=0, atol=2e-12)
    assert_allclose(actual_n.toarray(), reference_n.toarray(), rtol=0, atol=2e-12)
    assert len(actual_t) == len(reference_t)
    for actual, reference in zip(actual_t, reference_t):
        assert_allclose(actual.toarray(), reference.toarray(), rtol=0, atol=2e-12)
