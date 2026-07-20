#!/usr/bin/env python3
"""
NiO Anderson-impurity-model L-edge XAS using the experimental EDRIXS
SciPy pathway:

    setup_siam(...) -> ops(..., backend="scipy")
        -> ed_krylov_scipy(...)
        -> Lanczos continued-fraction XAS

This reproduces the physical model and output grid of:
    examples/sphinx/example_03_AIM_XAS.py

"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import numpy as np  # noqa: E402

import edrixs  # noqa: E402

from edrixs.solvers import (  # noqa: E402
    ed_krylov_scipy,
    ops,
    setup_siam,
    xas_krylov_scipy,
)


def build_physical_parameters() -> dict:
    """Build the NiO AIM parameters from the repository example."""
    nd = 8
    norb_d = 10
    nbath = 1
    v_noccu = nd + nbath * norb_d
    shell_name = ("d", "p")

    info = edrixs.utils.get_atom_data("Ni", "3d", nd, edge="L3")

    scale_dd = 0.8
    f2_dd = info["slater_i"][1][1] * scale_dd
    f4_dd = info["slater_i"][2][1] * scale_dd
    u_dd = 7.3
    f0_dd = u_dd + edrixs.get_F0("d", f2_dd, f4_dd)

    scale_dp = 0.8
    f2_dp = info["slater_n"][4][1] * scale_dp
    g1_dp = info["slater_n"][5][1] * scale_dp
    g3_dp = info["slater_n"][6][1] * scale_dp
    u_dp = 8.5
    f0_dp = u_dp + edrixs.get_F0("dp", g1_dp, g3_dp)

    slater = (
        [f0_dd, f2_dd, f4_dd],
        [f0_dd, f2_dd, f4_dd, f0_dp, f2_dp, g1_dp, g3_dp],
    )

    delta = 4.7
    e_d, e_l = edrixs.CT_imp_bath(u_dd, delta, nd)
    e_dc, e_lc, e_p = edrixs.CT_imp_bath_core_hole(
        u_dd,
        u_dp,
        delta,
        nd,
    )

    print(f"E_d  = {e_d:.6f} eV")
    print(f"E_L  = {e_l:.6f} eV")
    print(f"E_dc = {e_dc:.6f} eV")
    print(f"E_Lc = {e_lc:.6f} eV")
    print(f"E_p  = {e_p:.6f} eV")

    trans_c2n = edrixs.tmat_c2r("d", True)

    ten_dq = 0.56
    orbital_energies = np.repeat(
        np.array(
            [
                +0.6 * ten_dq,
                -0.4 * ten_dq,
                -0.4 * ten_dq,
                +0.6 * ten_dq,
                -0.4 * ten_dq,
            ]
        ),
        2,
    )
    crystal_field = np.diag(orbital_energies.astype(complex))

    # Preserve the exact published example: both impurity matrices use
    # the initial-state valence SOC value.
    zeta_d_i = info["v_soc_i"][0]
    soc_real = edrixs.cb_op(
        edrixs.atom_hsoc("d", zeta_d_i),
        trans_c2n,
    )
    imp_mat = crystal_field + soc_real + e_d * np.eye(norb_d)
    imp_mat_n = crystal_field + soc_real + e_dc * np.eye(norb_d)

    ten_dq_bath = 1.44
    bath_level = np.full((nbath, norb_d), e_l, dtype=complex)
    bath_level[0, :2] += ten_dq_bath * 0.6
    bath_level[0, 2:6] -= ten_dq_bath * 0.4
    bath_level[0, 6:8] += ten_dq_bath * 0.6
    bath_level[0, 8:] -= ten_dq_bath * 0.4

    bath_level_n = np.full((nbath, norb_d), e_lc, dtype=complex)
    bath_level_n[0, :2] += ten_dq_bath * 0.6
    bath_level_n[0, 2:6] -= ten_dq_bath * 0.4
    bath_level_n[0, 6:8] += ten_dq_bath * 0.6
    bath_level_n[0, 8:] -= ten_dq_bath * 0.4

    veg = 2.06
    vt2g = 1.21
    hyb = np.zeros((nbath, norb_d), dtype=complex)
    hyb[0, :2] = veg
    hyb[0, 2:6] = vt2g
    hyb[0, 6:8] = veg
    hyb[0, 8:] = vt2g

    om_shift = 857.6
    c_level = -om_shift - 5.0 * e_p

    return {
        "shell_name": shell_name,
        "nbath": nbath,
        "v_noccu": v_noccu,
        "slater": slater,
        "trans_c2n": trans_c2n,
        "imp_mat": imp_mat,
        "imp_mat_n": imp_mat_n,
        "bath_level": bath_level,
        "bath_level_n": bath_level_n,
        "hyb": hyb,
        "c_level": c_level,
        "c_soc": info["c_soc"],
        "om_shift": om_shift,
    }


def run(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    parameters = build_physical_parameters()

    temperature = 300.0
    neval = 50
    num_gs = 3
    thin = 0.0
    phi = 0.0

    ominc_xas = parameters["om_shift"] + np.linspace(-15.0, 25.0, 1000)
    gamma_c = np.full(ominc_xas.shape, 0.48 / 2.0)
    poltype_xas = [("isotropic", 0.0)]

    problem = setup_siam(
        parameters["shell_name"],
        parameters["nbath"],
        siam_type=0,
        v_noccu=parameters["v_noccu"],
        c_level=parameters["c_level"],
        c_soc=parameters["c_soc"],
        trans_c2n=parameters["trans_c2n"],
        imp_mat=parameters["imp_mat"],
        imp_mat_n=parameters["imp_mat_n"],
        bath_level=parameters["bath_level"],
        bath_level_n=parameters["bath_level_n"],
        hyb=parameters["hyb"],
        slater=parameters["slater"],
        ext_B=np.array([0.0, 0.0, 0.12]),
        on_which="spin",
        sparse_U=True,
    )

    hmat_i, hmat_n, trans_ops = ops(*problem, backend="scipy")

    print(f"Initial-space dimension:      {hmat_i.shape[0]}")
    print(f"Intermediate-space dimension: {hmat_n.shape[0]}")

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
        np.column_stack((np.arange(neval), eval_all, residuals)),
        header="state  energy_eV  residual_norm",
    )

    eval_i = eval_all[:num_gs]
    evec_i = evec_all[:, :num_gs]

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
        nkryl=200,
    )

    np.savetxt(
        output_dir / "xas.dat",
        np.column_stack((ominc_xas, xas)),
        header="incident_energy_eV  isotropic_xas",
    )

    fig, ax = plt.subplots()
    ax.plot(ominc_xas, xas[:, 0])
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("XAS intensity")
    ax.set_title("Anderson impurity model for NiO")
    fig.tight_layout()
    fig.savefig(output_dir / "xas.pdf")
    plt.close(fig)

    np.savez_compressed(
        output_dir / "nio_aim_xas_scipy_results.npz",
        eval_i=eval_all,
        residuals=residuals,
        ominc_xas=ominc_xas,
        gamma_c=gamma_c,
        xas=xas,
    )

    print(f"Results written to: {output_dir.resolve()}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the NiO AIM XAS benchmark through the SciPy pathway."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("NiO_AIM_XAS_scipy"),
        help="Directory for XAS, eigenvalues, plot, and NPZ results.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    run(arguments.output_dir)
