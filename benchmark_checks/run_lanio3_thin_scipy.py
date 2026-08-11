#!/usr/bin/env python3
"""LaNiO3 thin-film Ni L2,3-edge XAS/RIXS.

This example uses the EDRIXS backend-neutral solver interface with the
SciPy backend:

    setup_1v1c(...) -> get_ops(..., backend="scipy")
        -> ed(...) -> xas(...) / rixs(...)

This reproduces the physical model and output grids of:
    examples/more/RIXS/LaNiO3_thin/run_rixs_fsolver.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np  # noqa: E402


import edrixs  # noqa: E402

from edrixs.models import setup_1v1c  # noqa: E402
from edrixs.solvers import (  # noqa: E402
    get_ops,
    ed as solve_ed,
    xas as solve_xas,
    rixs as solve_rixs,
)


def build_physical_parameters():
    """Return the parameters used by the repository LaNiO3_thin example."""
    atom = edrixs.get_atom_data(
        "Ni",
        v_name="3d",
        v_noccu=8,
        edge="L23",
    )

    _, slater_i_atomic = [list(values) for values in zip(*atom["slater_i"])]
    _, slater_n_atomic = [list(values) for values in zip(*atom["slater_n"])]

    slater_i = edrixs.rescale(
        slater_i_atomic,
        ([1, 2], [0.65, 0.65]),
    )
    slater_i[0] = edrixs.get_F0("d", slater_i[1], slater_i[2])

    slater_n = edrixs.rescale(
        slater_n_atomic,
        ([1, 2, 4, 5, 6], [0.65, 0.65, 0.95, 0.7, 0.7]),
    )
    slater_n[0] = edrixs.get_F0("d", slater_n[1], slater_n[2])
    slater_n[3] = edrixs.get_F0("dp", slater_n[5], slater_n[6])

    zeta_d_i = atom["v_soc_i"][0]
    zeta_d_n = atom["v_soc_n"][0]
    zeta_p_n = (atom["edge_ene"][0] - atom["edge_ene"][1]) / 1.5

    crystal_field = edrixs.cf_tetragonal_d(
        1.3,   # 10Dq
        0.05,  # d1
        0.2,   # d3
    )

    return {
        "atom": atom,
        "slater": (slater_i, slater_n),
        "v_soc": (zeta_d_i, zeta_d_n),
        "c_soc": zeta_p_n,
        "crystal_field": crystal_field,
    }


def run(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    parameters = build_physical_parameters()
    atom = parameters["atom"]

    shell_name = ("d", "p")
    core_offset = 857.4
    temperature = 300.0
    neval = 10
    num_gs = 3

    thin = np.deg2rad(15.0)
    thout = np.deg2rad(75.0)
    phi = 0.0

    ominc_xas = np.linspace(core_offset - 10.0, core_offset + 20.0, 1000)
    ominc_rixs = np.linspace(core_offset - 5.9, core_offset - 0.9, 10)
    eloss = np.linspace(-0.5, 5.0, 1000)

    # L3-edge settings from the repository example.
    gamma_c = atom["gamma_c"][1]
    gamma_f = 0.1

    poltype_xas = [
        ("linear", 0.0),
        ("linear", np.pi / 2.0),
        ("left", 0.0),
        ("right", 0.0),
        ("isotropic", 0.0),
    ]
    poltype_rixs = [
        ("linear", 0.0, "linear", 0.0),
        ("linear", 0.0, "linear", np.pi / 2.0),
        ("linear", np.pi / 2.0, "linear", 0.0),
        ("linear", np.pi / 2.0, "linear", np.pi / 2.0),
    ]

    # -------------------------------------------------------------------------
    # 1. Define the orbital-space physical problem.
    # -------------------------------------------------------------------------
    problem = setup_1v1c(
        shell_name,
        shell_level=(0.0, -core_offset),
        v_soc=parameters["v_soc"],
        c_soc=parameters["c_soc"],
        v_noccu=8,
        slater=parameters["slater"],
        v_cfmat=parameters["crystal_field"],
        sparse_U=True,
    )

    # -------------------------------------------------------------------------
    # 2. Lift the problem into Fock space using SciPy operators.
    # -------------------------------------------------------------------------
    hmat_i, hmat_n, trans_ops = get_ops(
        *problem,
        backend="scipy",
    )

    print(f"Initial-space dimension:     {hmat_i.shape[0]}")
    print(f"Intermediate-space dimension: {hmat_n.shape[0]}")

    # -------------------------------------------------------------------------
    # 3. Calculate the Fortran eigenvalue count, then retain three states.
    # -------------------------------------------------------------------------
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
            "maxiter": 2000,
            "initial_guess": initial_guess,
            "suppress_lobpcg_warnings": False,
        },
    )

    residuals = np.array([
        np.linalg.norm(hmat_i @ evec_all[:, index] - eval_all[index] * evec_all[:, index])
        for index in range(neval)
    ])
    print("Computed initial energies:", eval_all)
    print("Eigenpair residual norms: ", residuals)

    np.savetxt(
        output_dir / "eval_i.dat",
        np.column_stack((np.arange(neval), eval_all, residuals)),
        header="state  energy_eV  residual_norm",
    )

    eval_i = eval_all[:num_gs]
    evec_i = evec_all[:, :num_gs]

    # The original Fortran input requests nkryl=100. The effective Krylov
    # dimensions cannot exceed the Hilbert-space dimensions (60 for XAS and
    # 45 for the final-state RIXS Lanczos chain).
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
            "nkryl": 100,
        },
    )

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
            "nkryl": min(100, hmat_i.shape[0]),
            "linsys_tol": 1.0e-11,
            "linsys_maxiter": 2000,
            "linsys_restart": min(60, hmat_n.shape[0]),
        },
    )

    # Match the files produced by the repository Fortran example.
    np.savetxt(
        output_dir / "xas.dat",
        np.column_stack((ominc_xas, xas)),
        header=(
            "energy_eV  linear_pi  linear_sigma  left_circular "
            "right_circular  isotropic"
        ),
    )

    rixs_pi = np.sum(rixs[:, :, 0:2], axis=2)
    rixs_sigma = np.sum(rixs[:, :, 2:4], axis=2)

    np.savetxt(
        output_dir / "rixs_pi.dat",
        np.column_stack((eloss, rixs_pi.T)),
        header="energy_loss_eV followed by one column per incident energy",
    )
    np.savetxt(
        output_dir / "rixs_sigma.dat",
        np.column_stack((eloss, rixs_sigma.T)),
        header="energy_loss_eV followed by one column per incident energy",
    )

    edrixs.plot_rixs_map(
        rixs_pi,
        ominc_rixs,
        eloss,
        str(output_dir / "rixsmap_pi.pdf"),
    )
    edrixs.plot_rixs_map(
        rixs_sigma,
        ominc_rixs,
        eloss,
        str(output_dir / "rixsmap_sigma.pdf"),
    )

    # A compact machine-readable bundle is useful for backend comparisons.
    np.savez_compressed(
        output_dir / "lanio3_thin_scipy_results.npz",
        eval_i=eval_all,
        residuals=residuals,
        ominc_xas=ominc_xas,
        xas=xas,
        ominc_rixs=ominc_rixs,
        eloss=eloss,
        rixs=rixs,
        rixs_pi=rixs_pi,
        rixs_sigma=rixs_sigma,
    )

    print(f"Results written to: {output_dir.resolve()}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the LaNiO3_thin benchmark through the backend-neutral "
            "interface using the SciPy backend."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("LaNiO3_thin_scipy"),
        help="Directory for spectra, maps, eigenvalues, and the NPZ result bundle.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    run(arguments.output_dir)
