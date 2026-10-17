"""Labeled factories for tests that start with explicitly constructed matrices."""

from functools import partial


def _identity(matrix):
    return matrix


def operator_factory(matrix, backend='scipy'):
    """Keep existing factory labels, or label an explicitly supplied matrix."""
    if hasattr(matrix, 'backend'):
        return matrix
    factory = partial(_identity, matrix)
    factory.backend = backend
    return factory


def operator_factories(matrices, backend='scipy'):
    return [operator_factory(matrix, backend) for matrix in matrices]


def materialize_ops(operators):
    """Obtain concrete matrices for numerical reference calculations."""
    hi, hn, transitions = operators
    if hi.backend == 'fortran':
        return operators
    return hi(), hn(), [transition() for transition in transitions]
