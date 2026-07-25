"""Consistency check for serial and process-based SciPy RIXS execution."""

import multiprocessing as mp

import pytest
from numpy.testing import assert_allclose

from edrixs.solvers import rixs_krylov_scipy

from ._helpers import exact_1v1c_reference_data

pytestmark = pytest.mark.integration


@pytest.mark.parallel
@pytest.mark.skipif(
    "fork" not in mp.get_all_start_methods(),
    reason="stable process-pool regression uses POSIX fork",
)
def test_parallel_rixs_matches_serial(small_1v1c_problem):
    """Compare final spectra from serial and two-process RIXS command chains.

    Both runs perform the same incoming transitions, GMRES solves, outgoing
    transitions, final-state responses, and broadening.  The only difference is
    whether independent contributions are evaluated serially or by two workers.
    """
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
        nkryl=hmat_i.shape[0],
        linsys_tol=1e-11,
        linsys_maxiter=200,
        linsys_restart=max(5, hmat_n.shape[0]),
    )

    serial = rixs_krylov_scipy(**common, workers=1)
    parallel = rixs_krylov_scipy(
        **common,
        workers=2,
        blas_threads=1,
        mp_start_method="fork",
    )

    assert_allclose(parallel, serial, rtol=1e-11, atol=1e-12)
