"""Option registry for the distributed PETSc/SLEPc backend."""

from __future__ import annotations

from numbers import Integral

from .._backend_options import (
    DynamicDefault,
    OptionSpec,
    boolean,
    nonnegative_real,
    optional_positive_integer,
    positive_integer,
    positive_real,
    validate_backend_kws,
)


def _communicator(value, label):
    """Accept PETSc communicators and mpi4py communicators."""
    protocols = (
        ('getRank', 'getSize'),
        ('Get_rank', 'Get_size'),
    )
    if not any(all(hasattr(value, method) for method in protocol)
               for protocol in protocols):
        raise TypeError(
            "{} must be a PETSc or mpi4py communicator".format(label)
        )


def _optional_nonnegative_integer(value, label):
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError("{} must be a nonnegative integer or None".format(label))
    if value < 0:
        raise ValueError("{} must be a nonnegative integer or None".format(label))


def _optional_string(value, label):
    if value is not None and not isinstance(value, str):
        raise TypeError("{} must be a string or None".format(label))


def _nonempty_string(value, label):
    if not isinstance(value, str):
        raise TypeError("{} must be a string".format(label))
    if not value:
        raise ValueError("{} must not be empty".format(label))


COMM = OptionSpec(
    DynamicDefault('PETSc.COMM_WORLD'),
    'PETSc or mpi4py communicator used to distribute the operator.',
    _communicator,
)

CONSTRUCTION_OPTIONS = {
    'comm': COMM,
    'tol_e': OptionSpec(
        1e-10,
        'Discard one-body coefficients whose magnitude does not exceed this threshold.',
        nonnegative_real,
    ),
    'tol_u': OptionSpec(
        1e-10,
        'Discard two-body coefficients whose magnitude does not exceed this threshold.',
        nonnegative_real,
    ),
    'nnz_guess_per_row': OptionSpec(
        None,
        'Nonzeros per row reserved separately for local and remote column blocks '
        'under MPI; estimated automatically when omitted.',
        _optional_nonnegative_integer,
    ),
    'mat_type': OptionSpec(
        None,
        "PETSc matrix type selected after assembly, for example 'aijcusparse'.",
        _optional_string,
    ),
    'assembly_chunk_cols': OptionSpec(
        4096,
        'Number of owned basis columns generated before entries are inserted.',
        positive_integer,
    ),
}

ED_OPTIONS = {
    'eigval_tol': OptionSpec(
        1e-8,
        'SLEPc residual tolerance: relative when shift is zero, absolute otherwise.',
        positive_real,
    ),
    'maxiter': OptionSpec(
        1000, 'Maximum SLEPc eigensolver iterations.', positive_integer
    ),
    'ncv': OptionSpec(
        None, 'SLEPc Krylov subspace size.', optional_positive_integer
    ),
    'verbose': OptionSpec(
        False, 'Print SLEPc convergence diagnostics on rank zero.', boolean
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
        1e-10,
        'PETSc KSP absolute tolerance (atol); shifted_pminres defaults to rtol=0.',
        positive_real,
    ),
    'linsys_maxiter': OptionSpec(
        1000, 'Maximum PETSc KSP iterations.', positive_integer
    ),
    'pminres_library': OptionSpec(
        None, 'Optional native library path for shifted_pminres.', _optional_string
    ),
    'ksp_type': OptionSpec(
        'shifted_pminres',
        'Shifted PMINRES (PCNONE, rtol=0), or a PETSc KSP solver type.',
        _nonempty_string
    ),
}


OPTIONS = {
    'build_op': CONSTRUCTION_OPTIONS,
    'get_ops': CONSTRUCTION_OPTIONS,
    'ed': ED_OPTIONS,
    'xas': XAS_OPTIONS,
    'rixs': RIXS_OPTIONS,
}


def validate_options(operation, backend_kws):
    """Validate options for a PETSc-backed operation."""
    return validate_backend_kws('petsc', operation, backend_kws, OPTIONS)
