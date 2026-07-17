#!/usr/bin/env python3
"""
URu2Si2 U O4,5-edge XAS/RIXS using the experimental EDRIXS SciPy pathway:

    setup_1v1c(...) -> ops(..., backend="scipy")
        -> ed_krylov_scipy(...)
        -> XAS/RIXS Krylov solvers

This reproduces the physical model and output grids of:
    examples/more/RIXS/URu2Si2/run_rixs_fsolver.py

"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np  # noqa: E402

import edrixs  # noqa: E402

from edrixs.solvers import (  # noqa: E402
    setup_1v1c,
    ops,
    ed_krylov_scipy,
    xas_krylov_scipy,
    rixs_krylov_scipy,
)


def build_physical_parameters() -> dict:
    """Return the parameters from the repository URu2Si2 Fortran example."""
    # 5f-5f Coulomb interactions, scaled to 77%.
    f2_ff = 9.711 * 0.77
    f4_ff = 6.364 * 0.77
    f6_ff = 4.677 * 0.77

    # The average 5f Coulomb interaction is set to zero in the fixed-occupancy
    # atomic calculation, exactly as in the reference example.
    uf_average = 0.0
    f0_ff = uf_average + edrixs.get_F0("f", f2_ff, f4_ff, f6_ff)

    # 5f-5d core-hole interactions, scaled to 60%.
    f2_fd = 10.652 * 0.60
    f4_fd = 6.850 * 0.60
    g1_fd = 12.555 * 0.60
    g3_fd = 7.768 * 0.60
    g5_fd = 5.544 * 0.60

    ufd_average = 0.0
    f0_fd = ufd_average + edrixs.get_F0("fd", g1_fd, g3_fd, g5_fd)

    slater = (
        [f0_ff, f2_ff, f4_ff, f6_ff],
        [
            f0_ff,
            f2_ff,
            f4_ff,
            f6_ff,
            f0_fd,
            f2_fd,
            f4_fd,
            g1_fd,
            g3_fd,
            g5_fd,
        ],
    )

    return {
        "slater": slater,
        "v_soc": (
            0.261 * 0.90,  # 5f SOC without a core hole
            0.274 * 0.90,  # 5f SOC with a core hole
        ),
        "c_soc": 3.2,      # 5d core SOC
    }


def build_gamma_c(ominc: np.ndarray) -> np.ndarray:
    """
    Reproduce the incident-energy-dependent core-hole broadening in the
    reference Fortran example.
    """
    gamma_c = np.zeros(len(ominc), dtype=float)
    begin = 101.5
    end = 104.0

    inside = np.count_nonzero((ominc > begin) & (ominc < end))
    if inside == 0:
        raise RuntimeError(
            "The incident-energy grid contains no points between "
            f"{begin} and {end} eV, so the reference broadening cannot be built."
        )

    step = (6.5 - 1.05) / (2.0 * inside)
    previous_count = 0

    for index, energy in enumerate(ominc):
        if energy <= begin:
            gamma_c[index] = 1.05 / 2.0
            previous_count += 1
        elif energy < end:
            gamma_c[index] = (
                1.05 / 2.0
                + step * (index - previous_count)
            )
        else:
            gamma_c[index] = 6.5 / 2.0

    return gamma_c


def run(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    parameters = build_physical_parameters()

    shell_name = ("f", "d")
    valence_occupancy = 2
    core_offset = 99.7
    temperature = 300.0
    neval = 30
    num_gs = 9

    thin = np.deg2rad(30.0)
    thout = np.deg2rad(60.0)
    phi = 0.0

    # The Fortran example uses the same 10-point incident-energy grid for
    # both XAS and RIXS.
    ominc = np.linspace(core_offset - 5.0, core_offset + 20.0, 10)
    eloss = np.linspace(-0.2, 2.5, 1000)

    gamma_c = build_gamma_c(ominc)
    gamma_f = 0.1

    poltype_xas = [
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
        v_noccu=valence_occupancy,
        slater=parameters["slater"],
        sparse_U=True,
    )

    # -------------------------------------------------------------------------
    # 2. Lift the model into the initial and intermediate Fock spaces.
    # -------------------------------------------------------------------------
    hmat_i, hmat_n, trans_ops = ops(
        *problem,
        backend="scipy",
    )

    print(f"Initial-space dimension:      {hmat_i.shape[0]}")
    print(f"Intermediate-space dimension: {hmat_n.shape[0]}")

    # -------------------------------------------------------------------------
    # 3. Compute 30 eigenpairs, then retain nine states for spectra.
    # -------------------------------------------------------------------------
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
        tol=1.0e-12,
        maxiter=3000,
        initial_guess=initial_guess,
        suppress_lobpcg_warnings=False,
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
        np.column_stack((
            np.arange(neval),
            eval_all,
            residuals,
        )),
        header="state  energy_eV  residual_norm",
    )

    eval_i = eval_all[:num_gs]
    evec_i = evec_all[:, :num_gs]

    np.savetxt(
        output_dir / "gamma_c.dat",
        np.column_stack((ominc, gamma_c)),
        header="incident_energy_eV  gamma_c_eV",
    )

    # -------------------------------------------------------------------------
    # 4. XAS.
    # -------------------------------------------------------------------------
    xas = xas_krylov_scipy(
        eval_i,
        evec_i,
        hmat_n,
        trans_ops,
        ominc,
        gamma_c=gamma_c,
        thin=thin,
        phi=phi,
        pol_type=poltype_xas,
        temperature=temperature,
        nkryl=100,
    )

    np.savetxt(
        output_dir / "xas.dat",
        np.column_stack((ominc, xas)),
        header="incident_energy_eV  isotropic_xas",
    )

    # -------------------------------------------------------------------------
    # 5. RIXS.
    # -------------------------------------------------------------------------
    rixs = rixs_krylov_scipy(
        eval_i,
        evec_i,
        hmat_i,
        hmat_n,
        trans_ops,
        ominc,
        eloss,
        gamma_c=gamma_c,
        gamma_f=gamma_f,
        thin=thin,
        thout=thout,
        phi=phi,
        pol_type=poltype_rixs,
        temperature=temperature,
        nkryl=100,
        linsys_tol=1.0e-10,
        linsys_maxiter=50000,
        linsys_restart=200,
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
        ominc,
        eloss,
        str(output_dir / "rixsmap_pi.pdf"),
    )
    edrixs.plot_rixs_map(
        rixs_sigma,
        ominc,
        eloss,
        str(output_dir / "rixsmap_sigma.pdf"),
    )

    np.savez_compressed(
        output_dir / "uru2si2_scipy_results.npz",
        eval_i=eval_all,
        residuals=residuals,
        ominc=ominc,
        gamma_c=gamma_c,
        xas=xas,
        eloss=eloss,
        rixs=rixs,
        rixs_pi=rixs_pi,
        rixs_sigma=rixs_sigma,
    )

    print(f"Results written to: {output_dir.resolve()}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the URu2Si2 benchmark through the new SciPy pathway."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("URu2Si2_scipy"),
        help="Directory for spectra, maps, eigenvalues, and the NPZ bundle.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    run(arguments.output_dir)
