.. _fortran-backend-options:

****************
Fortran keywords
****************

:doc:`Back to the backend overview <../backend>`

About this backend
==================

This backend runs the native, MPI-parallel Fortran solvers described in the
`EDRIXS paper <https://arxiv.org/abs/1812.05735>`__.  It uses a staged,
file-backed workflow: Python writes the solver inputs, calls the native
executable, and reads the results.  The
`EDRIXS repository <https://github.com/EDRIXS/edrixs>`__ contains the Fortran
sources and build instructions.

All staged operations are collective on an
`mpi4py communicator <https://mpi4py.github.io/mpi4py/stable/html/reference/mpi4py.MPI.Comm.html>`__.
When ``comm`` is omitted, ``mpi4py.MPI.COMM_WORLD`` is used.  Every rank in
that communicator must enter the operation; see the
`mpi4py collective communication tutorial <https://mpi4py.github.io/mpi4py/stable/html/tutorial.html#collective-communication>`__
for the underlying MPI model.

For exact diagonalization, ``ed_solver=0`` selects dense LAPACK,
``ed_solver=1`` selects Lanczos, and ``ed_solver=2`` selects parallel
ARPACK/PARPACK.  The `ARPACK project documentation
<https://www.arpack.org/installation>`__ describes the serial and MPI-enabled
libraries used by the last mode.

These options belong only to the staged interface.  The legacy ``*_fort``
functions retain their existing explicit parameters and behavior.

``get_ops``
===========

The Fortran backend constructs the complete problem at once, so
``build_op(..., backend='fortran')`` is not supported.

.. backend-options:: fortran get_ops

To resume a staged calculation from native input files already present in the
current working directory, call :func:`~edrixs.solvers.get_ops_disk`.  It
returns ``hmat_i``, ``hmat_n``, and ``trans_ops`` handles without rewriting
the problem files.

``ed``
======

.. backend-options:: fortran ed

For later XAS or RIXS calculations, ``idump`` must remain true and ``nvector``
must be at least the later ``num_gs`` value.  ``nvector`` cannot exceed
``num_evals``.

``xas``
=======

.. backend-options:: fortran xas

``rixs``
========

.. backend-options:: fortran rixs
