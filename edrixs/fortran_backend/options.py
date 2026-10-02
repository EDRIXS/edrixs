"""Option registry for the staged disk-backed Fortran interface."""

from __future__ import annotations

from .._backend_options import (
    DynamicDefault,
    OptionSpec,
    boolean,
    communicator,
    integer_choice,
    nonnegative_real,
    positive_integer,
    positive_real,
    validate_backend_kws,
)


COMM = OptionSpec(
    DynamicDefault('MPI.COMM_WORLD'),
    'mpi4py communicator used collectively by the backend.',
    communicator,
)

CONSTRUCTION_OPTIONS = {
    'comm': COMM,
    'tol': OptionSpec(
        1e-12,
        'Discard Coulomb coefficients whose magnitude does not exceed this threshold.',
        nonnegative_real,
    ),
}

ED_OPTIONS = {
    'comm': COMM,
    'ed_solver': OptionSpec(
        1,
        'Native eigensolver: 0 for dense, 1 for Lanczos, or 2 for ARPACK.',
        integer_choice(0, 1, 2),
    ),
    'nvector': OptionSpec(
        DynamicDefault('num_evals'),
        'Number of eigenvectors written for subsequent spectroscopy stages.',
        positive_integer,
    ),
    'ncv': OptionSpec(
        DynamicDefault('max(3, num_evals + 2)'),
        'ARPACK Krylov subspace size.',
        positive_integer,
    ),
    'idump': OptionSpec(True, 'Write eigenvectors to native files.', boolean),
    'maxiter': OptionSpec(
        500, 'Maximum eigensolver iterations.', positive_integer
    ),
    'min_ndim': OptionSpec(
        1000,
        'Dimension below which the native solver uses dense diagonalization.',
        positive_integer,
    ),
    'eigval_tol': OptionSpec(
        1e-8, 'Eigenvalue convergence tolerance.', positive_real
    ),
}

XAS_OPTIONS = {
    'comm': COMM,
    'num_gs': OptionSpec(
        DynamicDefault('len(eval_i)'),
        'Number of saved initial states included in the spectrum.',
        positive_integer,
    ),
    'nkryl': OptionSpec(
        200, 'Maximum intermediate-state Krylov dimension.', positive_integer
    ),
}

RIXS_OPTIONS = {
    'comm': COMM,
    'num_gs': OptionSpec(
        DynamicDefault('len(eval_i)'),
        'Number of saved initial states included in the spectrum.',
        positive_integer,
    ),
    'nkryl': OptionSpec(
        200, 'Maximum final-state Krylov dimension.', positive_integer
    ),
    'linsys_max': OptionSpec(
        1000,
        'Maximum iterations for the intermediate-state linear solve.',
        positive_integer,
    ),
    'linsys_tol': OptionSpec(
        1e-10,
        'Convergence tolerance for the intermediate-state linear solve.',
        positive_real,
    ),
}

OPTIONS = {
    'get_ops': CONSTRUCTION_OPTIONS,
    'ed': ED_OPTIONS,
    'xas': XAS_OPTIONS,
    'rixs': RIXS_OPTIONS,
}


def validate_options(operation, backend_kws):
    """Validate options for a staged Fortran backend operation."""
    return validate_backend_kws('fortran', operation, backend_kws, OPTIONS)
