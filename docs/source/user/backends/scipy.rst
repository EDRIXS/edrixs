.. _scipy-backend-options:

**************
SciPy keywords
**************

:doc:`Back to the backend overview <../backend>`

About this backend
==================

The ``scipy`` backend builds operators as SciPy
`CSR matrices <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html>`__.
The solver stages also accept a NumPy array, a SciPy sparse matrix, or a
matrix-free
`LinearOperator <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LinearOperator.html>`__.
ED uses SciPy's
`LOBPCG eigensolver <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lobpcg.html>`__;
XAS and RIXS use EDRIXS's Lanczos routines, and the RIXS correction-vector
step uses SciPy's
`GMRES solver <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gmres.html>`__.

Options are specific to the stage where they are supplied.

``build_op`` and ``get_ops``
============================

.. backend-options:: scipy get_ops

Here ``tol`` controls operator truncation.  It is separate from the
eigensolver tolerance accepted by ``ed``.

``ed``
======

.. backend-options:: scipy ed

``xas``
=======

.. backend-options:: scipy xas

``rixs``
========

.. backend-options:: scipy rixs
