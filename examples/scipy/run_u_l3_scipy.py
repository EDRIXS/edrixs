#!/usr/bin/env python3
"""U L3-edge XAS/RIXS using the EDRIXS backend-neutral solver interface with the SciPy backend:

    setup_2v1c(...) -> ops(..., backend="scipy")
        -> ed(...) -> xas(...) / rixs(...)

This reproduces the physical model and output grids of:
    examples/more/RIXS/U_L3/run_rixs_fsolver.py

Physical model
--------------
- Valence shells: U 5f and 6d
- Core shell: 2p3/2
- Transition: 2p3/2 -> 6d
- Total valence occupancy: 2
"""

from __future__ import annotations

import argparse
import collections
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np  # noqa: E402

import edrixs  # noqa: E402
from edrixs.solvers import (  # noqa: E402
    setup_2v1c,
    ops,
    ed as solve_ed,
    xas as solve_xas,
    rixs as solve_rixs,
)


def build_physical_parameters() -> dict:
    """Construct parameters used by the repository U_L3 Fortran example."""
    valence_occupancy = 2

    atom = edrixs.get_atom_data(
        "U",
        v_name=("5f", "6d"),
        v_noccu=(valence_occupancy, 0),
        edge="L3",
        trans_to_which=2,
        label=("f", "d", "p"),
    )

    slater_i = collections.OrderedDict(atom["slater_i"])
    slater_n = collections.OrderedDict(atom["slater_n"])

    slater_i["F2_ff"] *= 0.77
    slater_i["F4_ff"] *= 0.77
    slater_i["F6_ff"] *= 0.77
    slater_i["F0_ff"] = edrixs.get_F0(
        "f",
        slater_i["F2_ff"],
        slater_i["F4_ff"],
        slater_i["F6_ff"],
    )

    slater_n["F2_ff"] *= 0.77
    slater_n["F4_ff"] *= 0.77
    slater_n["F6_ff"] *= 0.77
    slater_n["F0_ff"] = edrixs.get_F0(
        "f",
        slater_n["F2_ff"],
        slater_n["F4_ff"],
        slater_n["F6_ff"],
    )

    slater_n["F0_fd"] = edrixs.get_F0(
        "fd",
        slater_n["G1_fd"],
        slater_n["G3_fd"],
        slater_n["G5_fd"],
    )
    slater_n["F0_fp"] = edrixs.get_F0(
        "fp",
        slater_n["G2_fp"],
        slater_n["G4_fp"],
    )
    slater_n["F0_dp"] = edrixs.get_F0(
        "dp",
        slater_n["G1_dp"],
        slater_n["G3_dp"],
    )

    return {
        "atom": atom,
        "valence_occupancy": valence_occupancy,
        "slater": (
            list(slater_i.values()),
            list(slater_n.values()),
        ),
        "v1_soc": (
            atom["v_soc_i"][0],
            atom["v_soc_n"][0],
        ),
        "v2_soc": (
            atom["v_soc_i"][1],
            atom["v_soc_n"][1],
        ),
    }


def run(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    parameters = build_physical_parameters()
    atom = parameters["atom"]

    shell_name = ("f", "d", "p32")
    temperature = 300.0

    # Retain the complete ninefold ground-state manifold for XAS/RIXS.
    neval = 20
    num_gs = 9

    thin = np.deg2rad(45.0)
    thout = np.deg2rad(45.0)
    phi = 0.0

    gamma_c = atom["gamma_c"][0]
    gamma_f = 0.1

    ominc_xas = np.linspace(-10.0, 20.0, 1000)
    ominc_rixs = np.linspace(0.0, 10.0, 10)
    eloss = np.linspace(-0.2, 5.0, 1000)

    poltype_xas = [("isotropic", 0.0)]
    poltype_rixs = [
        ("linear", 0.0, "linear", 0.0),
        ("linear", 0.0, "linear", np.pi / 2.0),
    ]

    # 1. Define the orbital-space model. trans_to_which=2 is essential:
    #    the L3 transition is 2p3/2 -> 6d, not 2p3/2 -> 5f.
    problem = setup_2v1c(
        shell_name,
        shell_level=(0.0, 5.0, 0.0),
        v1_soc=parameters["v1_soc"],
        v2_soc=parameters["v2_soc"],
        v_tot_noccu=parameters["valence_occupancy"],
        slater=parameters["slater"],
        trans_to_which=2,
        sparse_U=True,
    )

    # 2. Build many-body operators.
    hmat_i, hmat_n, trans_ops = ops(*problem, backend="scipy")

    print(f"Initial-space dimension:      {hmat_i.shape[0]}")
    print(f"Intermediate-space dimension: {hmat_n.shape[0]}")

    # 3. Compute 20 eigenpairs, then retain the ninefold ground manifold.
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
            "maxiter": 3000,
            "initial_guess": initial_guess,
            "suppress_lobpcg_warnings": False,
        },
    )

    residuals = np.array([
        np.linalg.norm(
            hmat_i @ evec_all[:, state]
            - eval_all[state] * evec_all[:, state]
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

    # 4. XAS.
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

    np.savetxt(
        output_dir / "xas.dat",
        np.column_stack((ominc_xas, xas)),
        header="incident_energy_eV  isotropic_xas",
    )

    # 5. RIXS.
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
            "linsys_tol": 1.0e-9,
            "linsys_maxiter": 1000,
            "linsys_restart": 200,
        },
    )

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

    np.savez_compressed(
        output_dir / "u_l3_scipy_results.npz",
        eval_i=eval_all,
        residuals=residuals,
        ominc_xas=ominc_xas,
        xas=xas,
        ominc_rixs=ominc_rixs,
        eloss=eloss,
        rixs=rixs,
        rixs_pi=rixs_pi,
        gamma_c=gamma_c,
    )

    print(f"Results written to: {output_dir.resolve()}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the U L3-edge 2v1c benchmark through the backend-neutral "
            "interface using the SciPy backend."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("U_L3_scipy"),
        help="Directory for spectra, map, eigenvalues, and the NPZ bundle.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    run(arguments.output_dir)
