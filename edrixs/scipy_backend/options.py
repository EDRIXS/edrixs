"""Option registry for the iterative SciPy backend."""

from __future__ import annotations

from .._backend_options import (
    OptionSpec,
    boolean,
    nonnegative_real,
    optional_positive_integer,
    positive_integer,
    positive_real,
    validate_backend_kws,
)


CONSTRUCTION_OPTIONS = {
    'tol': OptionSpec(
        1e-10,
        'Discard operator coefficients whose magnitude does not exceed this threshold.',
        nonnegative_real,
    ),
}

ED_OPTIONS = {
    'blocksize': OptionSpec(
        None,
        'LOBPCG block size; defaults to the requested number of eigenpairs.',
        optional_positive_integer,
    ),
    'tol': OptionSpec(1e-10, 'LOBPCG convergence tolerance.', positive_real),
    'maxiter': OptionSpec(200, 'Maximum LOBPCG iterations.', positive_integer),
    'seed': OptionSpec(None, 'Seed used to generate the initial LOBPCG block.'),
    'initial_guess': OptionSpec(
        None,
        'Initial LOBPCG block with shape (dimension, blocksize).',
    ),
    'suppress_lobpcg_warnings': OptionSpec(
        True,
        'Suppress SciPy nonconvergence warnings; results are still returned.',
        boolean,
    ),
}

XAS_OPTIONS = {
    'nkryl': OptionSpec(
        200, 'Maximum intermediate-state Lanczos dimension.', positive_integer
    ),
}

RIXS_OPTIONS = {
    'nkryl': OptionSpec(
        200, 'Maximum final-state Lanczos dimension.', positive_integer
    ),
    'linsys_tol': OptionSpec(
        1e-9, 'GMRES convergence tolerance.', positive_real
    ),
    'linsys_maxiter': OptionSpec(
        50000, 'Maximum GMRES iterations.', positive_integer
    ),
    'linsys_restart': OptionSpec(
        200, 'GMRES restart interval.', positive_integer
    ),
}

OPTIONS = {
    'build_op': CONSTRUCTION_OPTIONS,
    'get_ops': CONSTRUCTION_OPTIONS,
    'ed': ED_OPTIONS,
    'xas': XAS_OPTIONS,
    'rixs': RIXS_OPTIONS,
}


def validate_options(operation, backend_kws, *, backend='scipy'):
    """Validate options for a SciPy-backed operation."""
    return validate_backend_kws(backend, operation, backend_kws, OPTIONS)
