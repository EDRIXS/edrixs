#!/usr/bin/env python3
"""
LaNiO3 3d+4p / 2p two-valence-shell benchmark using the EDRIXS SciPy path.

This is the SciPy/Krylov counterpart of:
    examples/more/RIXS/LaNiO3_thin/test_2v1c.py

The model retains the original physical parameters and spectral grids, but uses:
    setup_2v1c -> ops -> ed_krylov_scipy
                 -> xas_krylov_scipy
                 -> rixs_krylov_scipy

The RIXS solver creates one process-pool job for every
(incident energy, polarization, retained initial state) triple.  By default it
uses all CPUs available to the process, capped by the number of jobs.
"""

import argparse
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np  # noqa: E402

import edrixs  # noqa: E402
from edrixs.solvers import (  # noqa: E402
    setup_2v1c,
    ops,
    ed_krylov_scipy,
    xas_krylov_scipy,
    rixs_krylov_scipy,
)


def build_physical_parameters():
    """Return the Slater integrals and one-body parameters of test_2v1c.py."""
    f2_dd = 12.234 * 0.65
    f4_dd = 7.598 * 0.65
    f0_dd = edrixs.get_F0("d", f2_dd, f4_dd)

    g1_d4p = 5.787 * 0.70
    g3_d4p = 3.291 * 0.70
    f0_d4p = edrixs.get_F0("dp", g1_d4p, g3_d4p)
    f2_d4p = 7.721 * 0.95

    f0_4p4p = 0.0
    f2_4p4p = 0.0

    f2_d2p = 7.721 * 0.95
    g1_d2p = 5.787 * 0.70
    g3_d2p = 3.291 * 0.70
    f0_d2p = edrixs.get_F0("dp", g1_d2p, g3_d2p)

    f0_4p2p = 2.0
    f2_4p2p = 2.0
    g0_4p2p = 2.0
    g2_4p2p = 2.0

    slater_initial = [
        f0_dd,
        f2_dd,
        f4_dd,
        f0_d4p,
        f2_d4p,
        g1_d4p,
        g3_d4p,
        f0_4p4p,
        f2_4p4p,
    ]
    slater_intermediate = [
        f0_dd,
        f2_dd,
        f4_dd,
        f0_d4p,
        f2_d4p,
        g1_d4p,
        g3_d4p,
        f0_4p4p,
        f2_4p4p,
        f0_d2p,
        f2_d2p,
        g1_d2p,
        g3_d2p,
        f0_4p2p,
        f2_4p2p,
        g0_4p2p,
        g2_4p2p,
    ]

    return {
        "slater": (slater_initial, slater_intermediate),
        "v1_soc": (0.083, 0.102),
        "c_soc": 11.24,
        "occupancy": 8,
        "crystal_field": edrixs.cf_tetragonal_d(1.3, 0.05, 0.2),
        "core_offset": 857.4,
    }


def run(
    output_dir,
    workers=None,
    blas_threads=1,
    mp_start_method=None,
    skip_xas=False,
    skip_rixs=False,
    sparse_u=True,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    timings = {}
    total_start = time.perf_counter()

    stage_start = time.perf_counter()
    parameters = build_physical_parameters()
    timings["parameter_setup_s"] = time.perf_counter() - stage_start

    shell_name = ("d", "p", "p")
    neval = 20
    num_gs = 3
    temperature = 300.0

    core_offset = parameters["core_offset"]
    thin = np.deg2rad(15.0)
    thout = np.deg2rad(75.0)
    phi = 0.0

    ominc_xas = np.linspace(core_offset - 10.0, core_offset + 20.0, 1000)
    ominc_rixs = np.linspace(core_offset - 5.9, core_offset - 0.9, 10)
    eloss = np.linspace(-0.5, 5.0, 1000)

    gamma_c = 0.2
    gamma_f = 0.1

    poltype_xas = [("isotropic", 0.0)]
    poltype_rixs = [
        ("linear", 0.0, "linear", 0.0),
        ("linear", 0.0, "linear", np.pi / 2.0),
    ]

    stage_start = time.perf_counter()
    problem = setup_2v1c(
        shell_name,
        shell_level=(0.0, 3.0, -core_offset),
        v1_soc=parameters["v1_soc"],
        c_soc=parameters["c_soc"],
        v1_cfmat=parameters["crystal_field"],
        v_tot_noccu=parameters["occupancy"],
        slater=parameters["slater"],
        trans_to_which=1,
        sparse_U=sparse_u,
    )
    timings["setup_2v1c_s"] = time.perf_counter() - stage_start

    stage_start = time.perf_counter()
    hmat_i, hmat_n, trans_ops = ops(*problem, backend="scipy")
    timings["ops_s"] = time.perf_counter() - stage_start

    print("Initial-space dimension:     ", hmat_i.shape[0])
    print("Intermediate-space dimension:", hmat_n.shape[0])

    stage_start = time.perf_counter()
    blocksize = neval
    rng = np.random.default_rng(12345)
    initial_guess = (
        rng.standard_normal((hmat_i.shape[0], blocksize))
        + 1j * rng.standard_normal((hmat_i.shape[0], blocksize))
    )
    eval_all, evec_all = ed_krylov_scipy(
        hmat_i,
        num_gs=neval,
        blocksize=blocksize,
        tol=1.0e-10,
        maxiter=1000,
        initial_guess=initial_guess,
        suppress_lobpcg_warnings=False,
    )
    timings["ed_s"] = time.perf_counter() - stage_start

    residuals = np.asarray([
        np.linalg.norm(hmat_i @ evec_all[:, i] - eval_all[i] * evec_all[:, i])
        for i in range(neval)
    ])

    print("Computed initial energies:")
    print(eval_all)
    print("Eigenpair residual norms:")
    print(residuals)

    np.savetxt(
        output_dir / "eval_i.dat",
        np.column_stack((np.arange(neval), eval_all, residuals)),
        header="state  energy_eV  residual_norm",
    )

    eval_i = eval_all[:num_gs]
    evec_i = evec_all[:, :num_gs]
    np.save(output_dir / "evec_i.npy", evec_i)

    xas = None
    if not skip_xas:
        stage_start = time.perf_counter()
        xas = xas_krylov_scipy(
            eval_i,
            evec_i,
            hmat_n,
            trans_ops,
            ominc_xas,
            gamma_c=gamma_c,
            thin=thin,
            phi=phi,
            pol_type=poltype_xas,
            temperature=temperature,
            nkryl=100,
        )
        timings["xas_s"] = time.perf_counter() - stage_start

        np.savetxt(
            output_dir / "xas.dat",
            np.column_stack((ominc_xas, xas)),
            header="incident_energy_eV  isotropic_xas",
        )

    rixs = None
    rixs_pi = None
    if not skip_rixs:
        stage_start = time.perf_counter()
        rixs = rixs_krylov_scipy(
            eval_i,
            evec_i,
            hmat_i,
            hmat_n,
            trans_ops,
            ominc_rixs,
            eloss,
            gamma_c=gamma_c,
            gamma_f=gamma_f,
            thin=thin,
            thout=thout,
            phi=phi,
            pol_type=poltype_rixs,
            temperature=temperature,
            nkryl=100,
            linsys_tol=1.0e-8,
            linsys_maxiter=500,
            linsys_restart=200,
            workers=workers,
            blas_threads=blas_threads,
            mp_start_method=mp_start_method,
        )
        timings["rixs_s"] = time.perf_counter() - stage_start

        rixs_pi = np.sum(rixs[:, :, 0:2], axis=2)
        np.savetxt(
            output_dir / "rixs_pi.dat",
            np.column_stack((eloss, rixs_pi.T)),
            header="energy_loss_eV followed by one column per incident energy",
        )
        edrixs.plot_rixs_map(
            rixs_pi,
            ominc_rixs,
            eloss,
            str(output_dir / "rixsmap_pi.pdf"),
        )

    timings["total_s"] = time.perf_counter() - total_start

    arrays = {
        "eval_i": eval_all,
        "residuals": residuals,
        "ominc_xas": ominc_xas,
        "ominc_rixs": ominc_rixs,
        "eloss": eloss,
        "gamma_c": np.asarray(gamma_c),
        "gamma_f": np.asarray(gamma_f),
    }
    if xas is not None:
        arrays["xas"] = xas
    if rixs is not None:
        arrays["rixs"] = rixs
        arrays["rixs_pi"] = rixs_pi

    np.savez_compressed(
        output_dir / "lanio3_2v1c_scipy_results.npz",
        **arrays
    )

    with (output_dir / "timings.json").open("w", encoding="utf-8") as stream:
        json.dump(timings, stream, indent=2, sort_keys=True)

    print("\nStage timings:")
    for name, seconds in timings.items():
        print("  {:22s} {:12.6f} s".format(name, seconds))
    print("\nResults written to:", output_dir.resolve())


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run the LaNiO3 3d+4p/2p two-valence-shell benchmark through "
            "the SciPy/Krylov pathway."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("LaNiO3_2v1c_scipy"),
        help="Directory for spectra, map, eigenpairs, timings, and NPZ data.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help=(
            "RIXS worker processes. The default uses all CPUs available to "
            "the process, capped by the 60 independent RIXS jobs."
        ),
    )
    parser.add_argument(
        "--blas-threads",
        type=int,
        default=1,
        help="Native BLAS/OpenMP threads per RIXS worker; keep at 1 normally.",
    )
    parser.add_argument(
        "--mp-start-method",
        choices=("fork", "spawn", "forkserver"),
        default=None,
        help="Optional multiprocessing start method.",
    )
    parser.add_argument(
        "--skip-xas",
        action="store_true",
        help="Run setup, operator construction, ED, and RIXS only.",
    )
    parser.add_argument(
        "--skip-rixs",
        action="store_true",
        help="Run setup, operator construction, ED, and XAS only.",
    )
    parser.add_argument(
        "--dense-u",
        action="store_true",
        help="Keep the orbital Coulomb tensors dense instead of sparse-flattened.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    run(
        args.output_dir,
        workers=args.workers,
        blas_threads=args.blas_threads,
        mp_start_method=args.mp_start_method,
        skip_xas=args.skip_xas,
        skip_rixs=args.skip_rixs,
        sparse_u=not args.dense_u,
    )


if __name__ == "__main__":
    main()
