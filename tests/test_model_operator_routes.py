"""Integration tests from model metadata through all operator-construction routes."""

import importlib.util
from itertools import product

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

from edrixs.fock_basis import FockBasisSpec
from edrixs.models import model_1v1c, model_2v1c, model_siam
from edrixs.solvers import get_ops


HAS_NUMBA = importlib.util.find_spec("numba") is not None
HAS_PETSC = importlib.util.find_spec("petsc4py") is not None


# Full public construction matrix: backend x basis representation x JIT choice.
ALL_ROUTES = [
    pytest.param(
        backend,
        basis_method,
        use_numba,
        id=f"{backend}-{basis_method}-{'numba' if use_numba else 'vanilla'}",
    )
    for backend, basis_method, use_numba in product(
        ("scipy", "petsc"), ("combinadic", "explicit"), (False, True)
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


def _to_dense(operator):
    """Convert small SciPy/NumPy/PETSc operators to a common dense reference."""
    if sp.issparse(operator):
        return operator.toarray()
    if isinstance(operator, np.ndarray):
        return operator

    # PETSc matrices in these tests use COMM_SELF and are therefore sequential.
    from petsc4py import PETSc

    nrow, ncol = operator.getSize()
    rows = np.arange(nrow, dtype=PETSc.IntType)
    cols = np.arange(ncol, dtype=PETSc.IntType)
    return np.asarray(operator.getValues(rows, cols))


def _route_kwargs(backend):
    """Keep PETSc tests sequential and force multiple bounded assembly chunks."""
    if backend == "petsc":
        from petsc4py import PETSc

        return {"backend_kws": {"comm": PETSc.COMM_SELF, "assembly_chunk_cols": 2}}
    return {}


def _skip_unavailable(backend, use_numba):
    """Skip only optional routes whose explicitly requested dependency is absent."""
    if use_numba and not HAS_NUMBA:
        pytest.skip("numba is required for the requested JIT route")
    if backend == "petsc" and not HAS_PETSC:
        pytest.skip("petsc4py is required for the PETSc route")


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
@pytest.mark.parametrize("backend,basis_method,use_numba", ALL_ROUTES)
def test_model_to_get_ops_all_eight_routes_are_numerically_consistent(
    backend, basis_method, use_numba
):
    """Compare all eight public operator-construction routes numerically.

    The explicit, non-Numba SciPy route is the common reference. Every other
    SciPy/PETSc x combinadic/explicit x vanilla/Numba route must produce the
    same initial Hamiltonian, intermediate Hamiltonian, and transition
    operators from the same compact model metadata.
    """
    _skip_unavailable(backend, use_numba)
    problem = _small_1v1c_problem()

    reference_i, reference_n, reference_t = get_ops(
        *problem,
        backend="scipy",
        basis_method="explicit",
        use_numba=False,
    )
    actual_i, actual_n, actual_t = get_ops(
        *problem,
        backend=backend,
        basis_method=basis_method,
        use_numba=use_numba,
        **_route_kwargs(backend),
    )

    assert_allclose(_to_dense(actual_i), reference_i.toarray(), rtol=0, atol=2e-12)
    assert_allclose(_to_dense(actual_n), reference_n.toarray(), rtol=0, atol=2e-12)
    assert len(actual_t) == len(reference_t)
    for actual, reference in zip(actual_t, reference_t):
        assert_allclose(_to_dense(actual), reference.toarray(), rtol=0, atol=2e-12)
