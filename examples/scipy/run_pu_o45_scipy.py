#!/usr/bin/env python3
"""Pu 5f6 O4,5-edge XAS/RIXS.

This example uses the EDRIXS backend-neutral solver interface with the
SciPy backend.

Reference: examples/more/RIXS/Pu_O45/run_rixs_fsolver.py
Expected dimensions: initial 3,003; intermediate 34,320.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np  # noqa: E402

import edrixs  # noqa: E402

from edrixs.solvers import (  # noqa: E402
    setup_1v1c,
    ops,
    ed as solve_ed,
    xas as solve_xas,
    rixs as solve_rixs,
)


def build_physical_parameters() -> dict:
    """Build the exact parameters from Pu_O45/run_rixs_fsolver.py."""
    occupancy = 6
    atom = edrixs.get_atom_data(
        "Pu", v_name="5f", v_noccu=occupancy, edge="O45"
    )

    _, slater_i_atomic = [list(values) for values in zip(*atom["slater_i"])]
    _, slater_n_atomic = [list(values) for values in zip(*atom["slater_n"])]

    slater_i = edrixs.rescale(
        slater_i_atomic,
        ([1, 2, 3], [0.77, 0.77, 0.77]),
    )
    slater_i[0] = edrixs.get_F0(
        "f", slater_i[1], slater_i[2], slater_i[3]
    )

    # Preserve the reference example exactly. It supplies eight indices but
    # seven factors, and edrixs.rescale uses zip(), so index 9 (G5_fd) remains
    # unscaled while indices 5, 6, 7, and 8 are scaled by 0.6.
    slater_n = edrixs.rescale(
        slater_n_atomic,
        (
            [1, 2, 3, 5, 6, 7, 8, 9],
            [0.77, 0.77, 0.77, 0.6, 0.6, 0.6, 0.6],
        ),
    )
    slater_n[0] = edrixs.get_F0(
        "f", slater_n[1], slater_n[2], slater_n[3]
    )
    slater_n[4] = edrixs.get_F0(
        "fd", slater_n[7], slater_n[8], slater_n[9]
    )

    zeta_f_i = atom["v_soc_i"][0] * 0.9
    zeta_f_n = atom["v_soc_n"][0] * 0.9
    zeta_d_n = (atom["edge_ene"][0] - atom["edge_ene"][1]) / 2.5
    om_shift = (
        2.0 * atom["edge_ene"][0] + 3.0 * atom["edge_ene"][1]
    ) / 5.0

    return {
        "occupancy": occupancy,
        "slater": (slater_i, slater_n),
        "v_soc": (zeta_f_i, zeta_f_n),
        "c_soc": zeta_d_n,
        "om_shift": om_shift,
    }


def run(
    output_dir: Path,
    *,
    skip_xas: bool = False,
    skip_rixs: bool = False,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    timings = {}
    total_start = time.perf_counter()

    stage_start = time.perf_counter()
    parameters = build_physical_parameters()
    timings["parameter_setup_s"] = time.perf_counter() - stage_start

    shell_name = ("f", "d")
    temperature = 300.0
    neval = 30
    num_gs = 1

    thin = np.deg2rad(45.0)
    thout = np.deg2rad(45.0)
    phi = 0.0

    om_shift = parameters["om_shift"]
    ominc_xas = np.linspace(om_shift - 5.0, om_shift + 20.0, 1000)
    ominc_rixs = np.linspace(om_shift - 5.0, om_shift + 20.0, 10)
    eloss = np.linspace(-0.2, 5.0, 1000)

    gamma_c = 0.2
    gamma_f = 0.1
    poltype_xas = [("isotropic", 0.0)]
    poltype_rixs = [
        ("linear", 0.0, "linear", 0.0),
        ("linear", 0.0, "linear", np.pi / 2.0),
    ]

    stage_start = time.perf_counter()
    problem = setup_1v1c(
        shell_name,
        shell_level=(0.0, -om_shift),
        v_soc=parameters["v_soc"],
        c_soc=parameters["c_soc"],
        v_noccu=parameters["occupancy"],
        slater=parameters["slater"],
        sparse_U=True,
    )
    timings["setup_1v1c_s"] = time.perf_counter() - stage_start

    stage_start = time.perf_counter()
    hmat_i, hmat_n, trans_ops = ops(*problem, backend="scipy")
    timings["ops_s"] = time.perf_counter() - stage_start

    print(f"Initial-space dimension:      {hmat_i.shape[0]}")
    print(f"Intermediate-space dimension: {hmat_n.shape[0]}")

    stage_start = time.perf_counter()
    blocksize = neval
    rng = np.random.default_rng(12345)
    initial_guess = (
        rng.standard_normal((hmat_i.shape[0], blocksize))
        + 1j * rng.standard_normal((hmat_i.shape[0], blocksize))
    )
    eval_all, evec_all = solve_ed(
        hmat_i,
        num_evals=neval,
        backend="scipy",
        backend_kws={
            "blocksize": blocksize,
            "tol": 1.0e-12,
            "maxiter": 5000,
            "initial_guess": initial_guess,
            "suppress_lobpcg_warnings": False,
        },
    )
    timings["ed_s"] = time.perf_counter() - stage_start

    residuals = np.array([
        np.linalg.norm(
            hmat_i @ evec_all[:, state] - eval_all[state] * evec_all[:, state]
        )
        for state in range(neval)
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

    xas = None
    if not skip_xas:
        stage_start = time.perf_counter()
        xas = solve_xas(
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
            backend="scipy",
            backend_kws={
                "nkryl": 200,
            },
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
        rixs = solve_rixs(
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
            backend="scipy",
            backend_kws={
            "parallel": True,
                "nkryl": 200,
                "linsys_tol": 1.0e-10,
                "linsys_maxiter": 5000,
                "linsys_restart": 200,
            },
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
        "gamma_c": gamma_c,
        "gamma_f": gamma_f,
    }
    if xas is not None:
        arrays["xas"] = xas
    if rixs is not None:
        arrays["rixs"] = rixs
        arrays["rixs_pi"] = rixs_pi

    np.savez_compressed(output_dir / "pu_o45_scipy_results.npz", **arrays)
    with (output_dir / "timings.json").open("w", encoding="utf-8") as stream:
        json.dump(timings, stream, indent=2, sort_keys=True)

    print("\nStage timings:")
    for name, seconds in timings.items():
        print(f"  {name:22s} {seconds:12.6f} s")
    print(f"\nResults written to: {output_dir.resolve()}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the Pu O4,5-edge benchmark through the backend-neutral "
            "interface using the SciPy backend."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("Pu_O45_scipy"),
        help="Directory for spectra, map, timing, and NPZ data.",
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
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    run(
        arguments.output_dir,
        skip_xas=arguments.skip_xas,
        skip_rixs=arguments.skip_rixs,
    )
