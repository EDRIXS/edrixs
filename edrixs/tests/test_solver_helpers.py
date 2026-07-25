"""Unit tests for internal handoffs in the SciPy sparse spectroscopy path.

These tests isolate conversions and per-contribution operations used between
``setup_*``, ``ops``, and the public XAS/RIXS drivers.  Complete comparisons of
final spectra live in ``solver_consistency_checks``.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.sparse.linalg import aslinearoperator

import edrixs._solvers_helpers as helpers


def test_dense_to_sparse_u_preserves_flattening_convention():
    """Check the Coulomb-data handoff from ``setup_*`` to sparse ``ops``.

    A dense four-index interaction may be flattened before the SciPy backend
    builds a many-electron Hamiltonian.  This test ensures orbital indices keep
    the convention expected by the sparse operator constructor.
    """
    rng = np.random.default_rng(3)
    umat = rng.normal(size=(3, 3, 3, 3))
    umat[np.abs(umat) < 0.6] = 0

    actual = helpers._umat_dense_to_sparse(umat).toarray()

    assert_allclose(actual, umat.reshape(9, 9))


def test_sparse_siam_embedding_matches_dense_embedding():
    """Check SIAM interaction placement before sparse Hamiltonian construction.

    ``setup_siam`` embeds impurity-core Coulomb terms into a larger space that
    also contains bath orbitals.  This test verifies the sparse embedding used
    by ``ops(..., backend="scipy")`` matches the dense placement exactly.
    """
    rng = np.random.default_rng(8)
    compact = rng.normal(size=(4, 4, 4, 4))
    dense = helpers._embed_impurity_core_umat(
        compact, v_norb=2, c_norb=2, ntot_v=6
    )
    sparse = helpers._embed_impurity_core_umat_sparse(
        compact, v_norb=2, c_norb=2, ntot_v=6
    )

    assert_allclose(sparse.toarray(), dense.reshape(64, 64))


def test_apply_linear_combination_handles_complex_and_all_zero_coefficients():
    """Check polarization-weighted transition application in XAS and RIXS.

    The public spectral drivers combine Cartesian transition operators using
    complex photon-polarization coefficients.  This test checks that stage,
    including the zero-polarization case, before any Krylov solve is started.
    """
    operators = [
        aslinearoperator(np.eye(2)),
        aslinearoperator(np.diag([2.0, 3.0])),
    ]
    vector = np.array([1.0, -2.0])

    actual = helpers._apply_linear_combination(operators, [1j, -0.5], vector)
    expected = 1j * vector - 0.5 * (np.diag([2.0, 3.0]) @ vector)

    assert_allclose(actual, expected)
    assert_allclose(
        helpers._apply_linear_combination(operators, [0, 0], vector),
        0,
    )


def test_expand_broadening_accepts_scalar_or_exact_length_array():
    """Check broadening preparation before XAS/RIXS spectrum assembly.

    The public drivers accept one lifetime width for a whole grid or one width
    per requested energy.  This test verifies that normalization step and
    rejects arrays that cannot be aligned with the final spectrum axes.
    """
    assert_allclose(helpers._expand_broadening(0.2, 3, "gamma"), [0.2] * 3)
    assert_allclose(
        helpers._expand_broadening([0.1, 0.2], 2, "gamma"),
        [0.1, 0.2],
    )

    with pytest.raises(ValueError, match="shape"):
        helpers._expand_broadening([0.1], 2, "gamma")


def test_zero_rhs_short_circuits_rixs_contribution(monkeypatch):
    """Check the incoming-transition stage of ``rixs_krylov_scipy``.

    Each RIXS contribution starts by applying the incoming photon operator to
    one retained initial state.  If that vector is zero, the contribution must
    return zero weight without entering the intermediate GMRES solve.
    """

    def forbidden_gmres(*args, **kwargs):
        raise AssertionError("GMRES must not run for a zero RHS")

    monkeypatch.setattr(helpers, "_gmres_scipy_compat", forbidden_gmres)
    record = helpers._rixs_krylov_one_contribution_scipy(
        hmat_i=aslinearoperator(np.eye(2)),
        hmat_n=aslinearoperator(np.eye(3)),
        trans_op_H=[],
        polvec_f=np.zeros(0),
        eval_i=np.array([0.0]),
        evec_i=np.array([[1.0], [0.0]]),
        istate=0,
        omega=1.0,
        gamma_c=0.2,
        rhs=np.zeros(3),
        skip_gs=False,
        nkryl=2,
        linsys_tol=1e-10,
        linsys_maxiter=10,
        linsys_restart=5,
    )

    assert record["norm"] == 0
    assert record["npoles"] == 1
    assert_allclose(record["alpha"], [0.0])
    assert record["beta"].size == 0


def test_gmres_failure_is_reported_with_context(monkeypatch):
    """Check failure reporting in the intermediate-state RIXS solve.

    ``rixs_krylov_scipy`` solves one correction-vector system for each incident
    energy and retained state.  This test ensures a failed solve identifies the
    affected contribution instead of silently contaminating the final spectrum.
    """
    monkeypatch.setattr(
        helpers,
        "_gmres_scipy_compat",
        lambda *args, **kwargs: (np.zeros(2, dtype=complex), 7),
    )

    with pytest.raises(RuntimeError, match=r"istate=0.*omega=1.25.*info=7"):
        helpers._rixs_krylov_one_contribution_scipy(
            hmat_i=aslinearoperator(np.eye(2)),
            hmat_n=aslinearoperator(np.eye(2)),
            trans_op_H=[aslinearoperator(np.eye(2))] * 3,
            polvec_f=np.array([1.0, 0.0, 0.0]),
            eval_i=np.array([0.0]),
            evec_i=np.array([[1.0], [0.0]]),
            istate=0,
            omega=1.25,
            gamma_c=0.2,
            rhs=np.ones(2),
            skip_gs=False,
            nkryl=2,
            linsys_tol=1e-10,
            linsys_maxiter=10,
            linsys_restart=5,
        )


def test_skip_gs_projects_every_retained_initial_state(monkeypatch):
    """Check elastic-state removal before the final RIXS response stage.

    After the intermediate solve and outgoing photon transition,
    ``rixs_krylov_scipy`` can remove all retained initial states when
    ``skip_gs=True``.  This test checks that projection before final poles and
    the energy-loss spectrum are constructed.
    """
    monkeypatch.setattr(
        helpers,
        "_gmres_scipy_compat",
        lambda *args, **kwargs: (np.array([2.0, -3.0], dtype=complex), 0),
    )

    common = dict(
        hmat_i=np.eye(2),
        hmat_n=np.eye(2),
        trans_op_H=[np.eye(2), np.zeros((2, 2)), np.zeros((2, 2))],
        polvec_f=np.array([1.0, 0.0, 0.0]),
        eval_i=np.array([0.0, 0.4]),
        evec_i=np.eye(2),
        istate=0,
        omega=1.0,
        gamma_c=0.2,
        rhs=np.ones(2),
        nkryl=2,
        linsys_tol=1e-10,
        linsys_maxiter=10,
        linsys_restart=5,
    )

    unprojected = helpers._rixs_krylov_one_contribution_scipy(
        **common,
        skip_gs=False,
    )
    projected = helpers._rixs_krylov_one_contribution_scipy(
        **common,
        skip_gs=True,
    )

    assert unprojected["norm"] == pytest.approx(13.0)
    assert projected["norm"] == 0.0
    assert projected["npoles"] == 1
    assert_allclose(projected["alpha"], [0.0])
    assert projected["beta"].size == 0
