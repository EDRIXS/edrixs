.. _dense-backend-options:

**************
Dense keywords
**************

:doc:`Back to the backend overview <../backend>`

About this backend
==================

The ``dense`` backend constructs NumPy
`ndarray objects <https://numpy.org/doc/stable/reference/arrays.ndarray.html>`__.
It uses SciPy's
`dense Hermitian eigensolver <https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html>`__
for ED.  XAS explicitly sums over all intermediate eigenstates, while RIXS
explicitly sums over all intermediate and final eigenstates.  These are the
same exact spectral-decomposition formulas used by the legacy dense Python
workflow.

The dense backend is useful for small problems, reference calculations, and
pedagogy.  It stores every matrix element and diagonalizes the intermediate
Hamiltonian during each XAS or RIXS call; RIXS also diagonalizes the final
Hamiltonian.  The sparse and matrix-free SciPy representations avoid that
storage cost and use iterative solvers for larger calculations.

Options are specific to the stage where they are supplied.

``build_op`` and ``get_ops``
============================

.. backend-options:: dense get_ops

Here ``tol`` controls operator truncation.

``ed``
======

.. backend-options:: dense ed

``xas``
========

.. backend-options:: dense xas

``rixs``
========

.. backend-options:: dense rixs
