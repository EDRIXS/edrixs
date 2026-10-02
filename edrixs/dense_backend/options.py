"""Option registry for the exact dense backend."""

from __future__ import annotations

from .._backend_options import validate_backend_kws
from ..scipy_backend.options import CONSTRUCTION_OPTIONS


OPTIONS = {
    'build_op': CONSTRUCTION_OPTIONS,
    'get_ops': CONSTRUCTION_OPTIONS,
    'ed': {},
    'xas': {},
    'rixs': {},
}


def validate_options(operation, backend_kws):
    """Validate options for an exact dense operation."""
    return validate_backend_kws('dense', operation, backend_kws, OPTIONS)
