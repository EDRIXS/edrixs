#!/usr/bin/env python3
"""Run the six SciPy examples sequentially and compare against legacy references.

Place this script in the directory containing the six ``run_*_scipy.py`` files
and the ``example_reference_data`` directory.  Each example is run with the
same Python interpreter that launched this script, so an activated conda EDRIXS
installation is used automatically.

The comparisons are intentionally compact:

* XAS/RIXS: energy grids are checked directly, and spectral values are compared
  with one normalized integrated difference.
* Eigenvalues: reference ``eigvals.dat`` energies are compared with the energy
  column of the SciPy example's ``eval_i.dat``.  The residual-norm column in
  ``eval_i.dat`` has no legacy counterpart and is not compared.

Default tolerances can be overridden from the command line.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np


EXAMPLES = (
    ("LaNiO3_2v1c", "run_lanio3_2v1c_scipy.py"),
     ("LaNiO3_thin", "run_lanio3_thin_scipy.py"),
     ("NiO_AIM_XAS", "run_nio_aim_xas_scipy.py"),
     ("U_L3", "run_u_l3_scipy.py"),
     ("URu2Si2", "run_uru2si2_scipy.py"),
    ("Pu_O45", "run_pu_o45_scipy.py"),
)

DEFAULT_SPECTRAL_TOLERANCE = 1.0e-6
DEFAULT_EIGENVALUE_ATOL = 1.0e-7
DEFAULT_GRID_ATOL = 1.0e-10


def _load_table(path: Path) -> np.ndarray:
    """Load a whitespace-separated numeric table and always return 2-D data."""
    data = np.loadtxt(path)
    return np.atleast_2d(data)


def _worst_error(actual: np.ndarray, reference: np.ndarray) -> tuple[float, float, tuple[int, ...]]:
    """Return maximum absolute error, scaled relative error, and worst index."""
    abs_error = np.abs(actual - reference)
    if abs_error.size == 0:
        return 0.0, 0.0, ()

    worst_flat = int(np.argmax(abs_error))
    worst_index = np.unravel_index(worst_flat, abs_error.shape)
    max_abs = float(abs_error[worst_index])

    scale = np.maximum(np.abs(reference), np.finfo(float).tiny)
    rel_error = abs_error / scale
    max_rel = float(np.max(rel_error))
    return max_abs, max_rel, worst_index


def _compare_eigenvalues(
    output_path: Path,
    reference_path: Path,
    *,
    atol: float,
) -> tuple[bool, str]:
    """Compare SciPy energies in ``eval_i.dat`` with legacy ``eigvals.dat``."""
    actual = _load_table(output_path)
    reference = _load_table(reference_path)

    if actual.shape[1] < 2 or reference.shape[1] < 2:
        return False, "expected at least two columns in both eigenvalue files"

    actual_energy = actual[:, 1]
    reference_energy = reference[:, 1]
    if actual_energy.shape != reference_energy.shape:
        return False, (
            "different number of eigenvalues: "
            f"actual={actual_energy.size}, reference={reference_energy.size}"
        )

    close = np.isclose(actual_energy, reference_energy, rtol=0.0, atol=atol)
    max_abs, _, worst = _worst_error(actual_energy, reference_energy)
    passed = bool(np.all(close))
    detail = f"max |ΔE|={max_abs:.3e} eV"
    if not passed:
        idx = int(worst[0])
        detail += (
            f"; worst state={idx + 1}, actual={actual_energy[idx]:.12g}, "
            f"reference={reference_energy[idx]:.12g}"
        )
    return passed, detail


def _compare_spectrum(
    output_path: Path,
    reference_path: Path,
    *,
    tolerance: float,
    grid_atol: float,
) -> tuple[bool, str]:
    """Compare grids directly and spectra by normalized integrated difference."""
    actual = _load_table(output_path)
    reference = _load_table(reference_path)

    if actual.shape != reference.shape:
        return False, f"shape mismatch: actual={actual.shape}, reference={reference.shape}"
    if actual.shape[1] < 2:
        return False, "expected energy grid plus at least one spectrum column"

    actual_grid = actual[:, 0]
    reference_grid = reference[:, 0]
    grid_close = np.isclose(actual_grid, reference_grid, rtol=0.0, atol=grid_atol)
    grid_max_abs, _, grid_worst = _worst_error(actual_grid, reference_grid)
    if not np.all(grid_close):
        row = int(grid_worst[0])
        return False, (
            f"energy-grid mismatch: max |Δx|={grid_max_abs:.3e} at row {row + 1}; "
            f"actual={actual_grid[row]:.12g}, reference={reference_grid[row]:.12g}"
        )

    actual_values = actual[:, 1:]
    reference_values = reference[:, 1:]
    numerator = float(np.sum(np.abs(actual_values - reference_values)))
    denominator = float(np.sum(np.abs(actual_values + reference_values)))
    if denominator == 0.0:
        spectral_difference = 0.0 if numerator == 0.0 else np.inf
    else:
        spectral_difference = numerator / denominator

    passed = bool(spectral_difference <= tolerance)
    detail = (
        f"grid max |Δ|={grid_max_abs:.3e}; "
        f"normalized spectral difference={spectral_difference:.3e}"
    )
    return passed, detail


def _expected_comparisons(reference_dir: Path) -> list[tuple[str, Path, Path]]:
    """Return comparisons required by the files present in one reference case."""
    comparisons = [
        ("eigenvalues", Path("eval_i.dat"), reference_dir / "eigvals.dat"),
    ]
    spectrum_names = ["xas.dat", "rixs_pi.dat", "rixs_sigma.dat"]
    if reference_dir.name == "LaNiO3_2v1c":
        spectrum_names.remove("rixs_sigma.dat")
    for name in spectrum_names:
        reference = reference_dir / name
        if reference.exists():
            comparisons.append((name, Path(name), reference))
    return comparisons


def _run_example(script: Path, output_dir: Path) -> tuple[bool, float]:
    """Run one example with the current Python interpreter."""
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    command = [sys.executable, str(script), "--output-dir", str(output_dir)]

    print(f"\nRUN  {script.name}")
    print("     " + " ".join(command))
    started = time.perf_counter()
    result = subprocess.run(command, cwd=script.parent, env=env, check=False)
    elapsed = time.perf_counter() - started
    if result.returncode != 0:
        print(f"FAIL {script.name}: exit code {result.returncode} ({elapsed:.1f} s)")
        return False, elapsed

    print(f"DONE {script.name} ({elapsed:.1f} s)")
    return True, elapsed


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help=(
            "directory for fresh example outputs "
            "(default: example_validation_output beside this script)"
        ),
    )
    parser.add_argument(
        "--compare-only",
        action="store_true",
        help="do not run examples; compare files already present under --output-root",
    )
    parser.add_argument(
        "--spectral-tolerance",
        type=float,
        default=DEFAULT_SPECTRAL_TOLERANCE,
        help=(
            "maximum normalized integrated XAS/RIXS difference "
            f"(default: {DEFAULT_SPECTRAL_TOLERANCE:g})"
        ),
    )
    parser.add_argument(
        "--eigenvalue-atol",
        type=float,
        default=DEFAULT_EIGENVALUE_ATOL,
        help=f"absolute eigenvalue tolerance in eV (default: {DEFAULT_EIGENVALUE_ATOL:g})",
    )
    parser.add_argument(
        "--grid-atol",
        type=float,
        default=DEFAULT_GRID_ATOL,
        help=f"absolute tolerance for XAS/RIXS energy grids (default: {DEFAULT_GRID_ATOL:g})",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    base_dir = Path(__file__).resolve().parent
    reference_root = base_dir / "example_reference_data"
    output_root = args.output_root or (base_dir / "example_validation_output")
    if not output_root.is_absolute():
        output_root = (Path.cwd() / output_root).resolve()

    missing = []
    for case_name, script_name in EXAMPLES:
        script = base_dir / script_name
        reference_dir = reference_root / case_name
        if not script.exists():
            missing.append(str(script))
        if not reference_dir.is_dir():
            missing.append(str(reference_dir))
    if missing:
        print("Missing required paths:", file=sys.stderr)
        for path in missing:
            print(f"  {path}", file=sys.stderr)
        return 2

    print("SciPy example validation")
    print(f"Python:          {sys.executable}")
    print(f"Examples:        {base_dir}")
    print(f"References:      {reference_root}")
    print(f"Outputs:         {output_root}")
    print(
        "Tolerances:      "
        f"spectral difference={args.spectral_tolerance:g}; "
        f"eigenvalue atol={args.eigenvalue_atol:g} eV; grid atol={args.grid_atol:g}"
    )

    run_status: dict[str, bool] = {}
    elapsed_times: dict[str, float] = {}
    if not args.compare_only:
        output_root.mkdir(parents=True, exist_ok=True)
        for case_name, script_name in EXAMPLES:
            ok, elapsed = _run_example(base_dir / script_name, output_root / case_name)
            run_status[case_name] = ok
            elapsed_times[case_name] = elapsed
    else:
        run_status = {case_name: True for case_name, _ in EXAMPLES}

    print("\n" + "=" * 80)
    print("COMPARISON RESULTS")
    print("=" * 80)

    overall_pass = True
    for case_name, _ in EXAMPLES:
        print(f"\n{case_name}")
        if not run_status[case_name]:
            print("  FAIL  example execution failed; comparison skipped")
            overall_pass = False
            continue

        output_dir = output_root / case_name
        reference_dir = reference_root / case_name
        case_pass = True
        for label, output_name, reference_path in _expected_comparisons(reference_dir):
            output_path = output_dir / output_name
            if not output_path.exists():
                print(f"  FAIL  {label}: missing output file {output_path.name}")
                case_pass = False
                continue

            if label == "eigenvalues":
                passed, detail = _compare_eigenvalues(
                    output_path,
                    reference_path,
                    atol=args.eigenvalue_atol,
                )
            else:
                passed, detail = _compare_spectrum(
                    output_path,
                    reference_path,
                    tolerance=args.spectral_tolerance,
                    grid_atol=args.grid_atol,
                )

            status = "PASS" if passed else "FAIL"
            print(f"  {status}  {label}: {detail}")
            case_pass &= passed

        if not args.compare_only:
            print(f"  runtime: {elapsed_times[case_name]:.1f} s")
        overall_pass &= case_pass

    print("\n" + "=" * 80)
    if overall_pass:
        print("PASS: all six examples agree with the legacy reference data.")
        return 0

    print("FAIL: one or more examples differ from the legacy reference data.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
