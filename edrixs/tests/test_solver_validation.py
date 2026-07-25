"""Unit tests for public input contracts of SciPy XAS and RIXS.

These tests check the boundary after Hamiltonians and transition operators have
been prepared but before any response calculation starts.  They ensure shape,
broadening, polarization, solver, and multiprocessing errors are reported early.
"""

import numpy as np
import pytest
from scipy.sparse.linalg import LinearOperator

import edrixs._solvers_helpers as helpers
from edrixs.solvers import rixs_krylov_scipy, xas_krylov_scipy


def _valid_xas_kwargs():
    """Return the smallest valid public XAS call used as a validation baseline."""
    return dict(
        eval_i=np.array([0.0]),
        evec_i=np.array([[1.0], [0.0]], dtype=complex),
        hmat_n=np.eye(3),
        trans_op=[np.zeros((3, 2), dtype=complex) for _ in range(3)],
        ominc=np.array([0.8, 1.0]),
        gamma_c=0.2,
        nkryl=3,
    )


@pytest.mark.parametrize(
    "update, message",
    [
        ({"eval_i": np.array([[0.0]])}, "one-dimensional"),
        ({"evec_i": np.array([1.0, 0.0])}, "two-dimensional"),
        ({"evec_i": np.eye(2, dtype=complex)}, "columns"),
        ({"hmat_n": np.zeros((3, 2))}, "square"),
        ({"trans_op": [np.zeros((3, 2)) for _ in range(2)]}, r"len\(trans_op\)"),
        (
            {"trans_op": [np.zeros((3, 2)), np.zeros((3, 2)), np.zeros((2, 2))]},
            "expected",
        ),
        ({"scatter_axis": np.eye(2)}, "scatter_axis"),
        ({"nkryl": 0}, "positive"),
        ({"gamma_c": [0.2]}, "shape"),
    ],
)
def test_xas_public_validation(update, message):
    """Reject malformed data before ``xas_krylov_scipy`` builds a response.

    The public XAS stage receives retained initial states, an intermediate
    Hamiltonian, transitions, geometry, and an energy grid.  This test checks
    their shared contract before start vectors and Krylov poles are generated.
    """
    kwargs = _valid_xas_kwargs()
    kwargs.update(update)
    with pytest.raises(ValueError, match=message):
        xas_krylov_scipy(**kwargs)


def test_xas_rejects_unknown_polarization():
    """Reject an unsupported photon polarization before XAS transitions are applied.

    Polarization coefficients select and combine transition operators at the
    start of the XAS command chain.  This test prevents an undefined geometry
    from reaching the intermediate-state response and final spectrum.
    """
    kwargs = _valid_xas_kwargs()
    kwargs["pol_type"] = [("diagonal", 0.0)]
    with pytest.raises(ValueError, match="Unknown XAS polarization"):
        xas_krylov_scipy(**kwargs)


def _valid_rixs_kwargs():
    """Return the smallest valid public RIXS call used as a validation baseline."""
    transitions = [np.zeros((3, 2), dtype=complex) for _ in range(3)]
    transitions[0][0, 0] = 1.0
    return dict(
        eval_i=np.array([0.0]),
        evec_i=np.array([[1.0], [0.0]], dtype=complex),
        hmat_i=np.eye(2),
        hmat_n=np.eye(3),
        trans_op=transitions,
        ominc=np.array([1.0]),
        eloss=np.array([0.0, 0.2]),
        gamma_c=0.2,
        gamma_f=0.1,
        nkryl=2,
        linsys_maxiter=20,
        linsys_restart=5,
        workers=1,
    )


@pytest.mark.parametrize(
    "update, message",
    [
        ({"hmat_i": np.zeros((2, 3))}, "hmat_i must be square"),
        ({"hmat_n": np.zeros((3, 2))}, "hmat_n must be square"),
        ({"eval_i": np.array([[0.0]])}, "one-dimensional"),
        ({"evec_i": np.array([1.0, 0.0])}, "two-dimensional"),
        ({"eval_i": np.array([0.0, 0.2])}, r"len\(eval_i\) must equal"),
        ({"evec_i": np.ones((3, 1))}, r"evec_i.shape\[0\]"),
        ({"trans_op": [np.zeros((3, 2)) for _ in range(2)]}, r"len\(trans_op\)"),
        (
            {"trans_op": [np.zeros((3, 2)), np.zeros((3, 2)), np.zeros((2, 2))]},
            "expected",
        ),
        ({"scatter_axis": np.eye(2)}, "scatter_axis"),
        ({"nkryl": 0}, "nkryl"),
        ({"linsys_maxiter": 0}, "linsys_maxiter"),
        ({"linsys_restart": 0}, "linsys_restart"),
        ({"blas_threads": 0}, "blas_threads"),
        ({"workers": 0}, "workers"),
        ({"gamma_c": [0.1, 0.2]}, "gamma_c.*shape"),
        ({"gamma_f": [0.1]}, "gamma_f.*shape"),
    ],
)
def test_rixs_public_validation(update, message):
    """Reject malformed data before ``rixs_krylov_scipy`` starts its command chain.

    RIXS combines initial states, two Hamiltonians, transition operators, two
    energy grids, linear-solver controls, and worker settings.  This test checks
    those inputs before incoming transitions, GMRES solves, and final poles.
    """
    kwargs = _valid_rixs_kwargs()
    kwargs.update(update)
    with pytest.raises(ValueError, match=message):
        rixs_krylov_scipy(**kwargs)


def test_rixs_requires_transition_adjoint_action():
    """Require the outgoing-transition action used after the RIXS GMRES solve.

    RIXS first applies a transition into the intermediate space and later needs
    its adjoint to return to the final-state space.  This test rejects operators
    that cannot perform that second stage before any spectrum is calculated.
    """
    no_adjoint = LinearOperator(
        shape=(3, 2),
        matvec=lambda vector: np.zeros(3, dtype=complex),
        dtype=np.complex128,
    )
    kwargs = _valid_rixs_kwargs()
    kwargs["trans_op"] = [no_adjoint, no_adjoint, no_adjoint]

    with pytest.raises(TypeError, match="adjoint action"):
        rixs_krylov_scipy(**kwargs)


def test_process_pool_rejects_nonserializable_operator():
    """Check operator preparation for multicore ``rixs_krylov_scipy``.

    Before independent RIXS contributions are sent to worker processes, their
    Hamiltonians and transitions must be transferable.  This test reports an
    unsupported operator before workers start and before partial spectra exist.
    """
    nonserializable = LinearOperator(
        shape=(2, 2),
        matvec=lambda vector: vector,
        dtype=np.complex128,
    )

    with pytest.raises(TypeError, match="not serializable"):
        helpers._operator_for_process_pool(nonserializable, "custom_op")


def test_rixs_worker_job_requires_initializer(monkeypatch):
    """Check worker-state initialization in multicore RIXS execution.

    Each process receives Hamiltonians, transitions, states, and solver settings
    once before it evaluates incident-energy contributions.  This test ensures
    a worker cannot run without that state and return incomplete pole records.
    """
    monkeypatch.setattr(helpers, "_RIXS_WORKER_STATE", None)
    with pytest.raises(RuntimeError, match="not initialized"):
        helpers._rixs_pool_job((0, 0, 0))
