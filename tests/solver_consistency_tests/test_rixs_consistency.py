"""Consistency check for explicit and inferred SciPy RIXS dispatch."""

import pytest
from numpy.testing import assert_allclose

from edrixs.solvers import rixs

from ._helpers import exact_1v1c_reference_data


pytestmark = pytest.mark.integration


def test_explicit_and_inferred_scipy_rixs_match(small_1v1c_problem):
    """Compare explicit and inferred backend dispatch through SciPy RIXS."""
    hmat_i, hmat_n, transitions, eval_i, evec_i, _, _ = (
        exact_1v1c_reference_data(small_1v1c_problem)
    )
    common = dict(
        eval_i=eval_i[:1],
        evec_i=evec_i[:, :1],
        hmat_i=hmat_i,
        hmat_n=hmat_n,
        trans_op=transitions,
        ominc=[0.8, 1.0],
        eloss=[0.0, 0.2, 0.4],
        gamma_c=0.2,
        gamma_f=0.1,
        backend_kws={
            "nkryl": hmat_i.shape[0],
            "linsys_tol": 1e-11,
            "linsys_maxiter": 200,
            "linsys_restart": max(5, hmat_n.shape[0]),
        },
    )

    explicit = rixs(**common, backend="scipy")
    inferred = rixs(**common)

    assert_allclose(inferred, explicit, rtol=1e-11, atol=1e-12)
