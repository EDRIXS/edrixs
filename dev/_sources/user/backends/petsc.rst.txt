.. _petsc-backend-options:

**************
PETSc keywords
**************

:doc:`Back to the backend overview <../backend>`

About this backend
==================

PETSc supplies distributed vectors and matrices together with scalable linear
solvers and preconditioners; its concise
`architecture overview <https://petsc.org/release/overview/nutshell/>`__
introduces the ``Vec``, ``Mat``, ``KSP``, and ``PC`` objects.  Python programs
access them through the official
`petsc4py API <https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.html>`__.
For eigenvalue problems, SLEPc builds on PETSc and exposes its ``EPS`` solver
through the
`SLEPc eigenvalue guide <https://slepc.upv.es/release/documentation/manual/eps.html>`__
and the
`slepc4py EPS API <https://slepc.upv.es/release/slepc4py/reference/slepc4py.SLEPc.EPS.html>`__.

EDRIXS constructs distributed PETSc ``Mat`` objects, obtains low-energy
eigenpairs with SLEPc's Krylov--Schur method, and evaluates XAS and RIXS with
distributed Lanczos and KSP solves.  Importing EDRIXS does not import PETSc;
``petsc4py`` and ``slepc4py`` are loaded only when this backend is used.

``build_op`` and ``get_ops`` return factories labeled with
``.backend == 'petsc'``. Calling a factory assembles its distributed matrix;
public ``ed``, ``xas``, and ``rixs`` read the labels and pass the factories to
the PETSc backend. All MPI ranks in the matrix communicator must participate
in construction and solver calls. See :ref:`lazy-operators`.

Options are specific to the stage where they are supplied.

Matrix lifetime
===============

PETSc ``xas`` and ``rixs`` accept ``matrix_lifetime`` in ``backend_kws``:

* ``'retain'`` (default) builds each required factory on first use, keeps its
  matrix for reuse during the calculation, and destroys it before returning.
* ``'recompute'`` builds a matrix for each operation and destroys it immediately
  afterward, rebuilding it for later uses. This reduces peak memory usage at the
  cost of more CPU time. This can be useful for the largest jobs.


``build_op`` and ``get_ops``
============================

.. backend-options:: petsc get_ops

Two matrix construction options are available. ``assembly_method='row'`` is
the default and inserts entries by row in bounded chunks.
``assembly_method='coo'`` counts contributions in a first pass, then generates
them again into exact-sized coordinate and value arrays for bulk PETSc COO
assembly. PETSc combines duplicate coordinates and communicates remote entries.
Both options produce AIJ (CSR) matrices before optional ``mat_type`` conversion;
COO describes the assembly interface, rather than the final storage format.
Both options support the same solver lifetime policies.
COO avoids the Python row-insertion loop but requires full local entry buffers
and PETSc's COO mappings, so compare assembly time and peak memory for your
problem. ``assembly_chunk_cols`` bounds individual kernel calls for both
methods; COO still collects all local contributions before insertion.

``nnz_guess_per_row`` controls PETSc matrix preallocation.  Leaving it unset
uses an estimate from the retained one- and two-body coefficients. With multiple
MPI ranks, the hint is applied separately to the local-column and remote-column
blocks of each rank's matrix rows. With one rank it applies to the sequential
matrix. Larger hints reserve more memory but can avoid reallocations during
row assembly; they do not change the matrix values. COO computes its pattern
directly and ignores this hint. A
``mat_type`` such as ``'aijcusparse'`` can select an accelerated matrix format;
see PETSc's `matrix type overview
<https://petsc.org/release/manual/mat/#basic-matrix-operations>`__.

``ed``
======

Pass the final model output as ``ed(hmat_i, shift=shift)`` to subtract the
filled-core shell energy from the factory-created distributed Hamiltonian before
diagonalization. Returned eigenvalues retain the original energy reference.
With a nonzero ``shift``, ``eigval_tol`` is an absolute residual tolerance for
the shifted problem, avoiding division by a near-zero centered eigenvalue.
With ``shift=0``, it is a relative residual tolerance, scaled by the eigenvalue
magnitude.

ED calls its factory once and destroys the matrix before returning.
``matrix_lifetime`` is not accepted for ED.

.. backend-options:: petsc ed

The eigensolver uses SLEPc ``EPS`` with a Hermitian problem type and requests
the smallest-real eigenpairs.  ``ncv`` sets the Krylov subspace size described
by `EPS.setDimensions
<https://slepc.upv.es/release/slepc4py/reference/slepc4py.SLEPc.EPS.html#slepc4py.SLEPc.EPS.setDimensions>`__.

``xas``
=======

XAS accepts ``matrix_lifetime`` in ``backend_kws``. The default ``'retain'``
builds required transition components and the intermediate Hamiltonian on first
use and keeps them until calculation exit. With ``'recompute'``, each transition
component is destroyed after its application, and the intermediate Hamiltonian
is rebuilt and destroyed for each nonzero seed's Lanczos calculation. At most
one factory-created matrix is alive at a time in this mode. Components with
zero polarization coefficients and Hamiltonians for zero seeds are skipped.
Both policies destroy factory-created matrices before returning, including on
failure; initial-state eigenvectors remain owned by the caller.

XAS and RIXS center factory-created distributed matrices using the lowest retained initial
energy and the first incident energy, which should lie near resonance. RIXS
also centers the final Hamiltonian. Spectra are evaluated on these relative
energies before restoring the original reference in returned pole dictionaries.
Factory-created Hamiltonians are shifted in place.

.. backend-options:: petsc xas

``rixs``
========

.. backend-options:: petsc rixs

The default solver is ``shifted_pminres``. ``linsys_tol`` sets KSP's absolute
tolerance, ``atol``, and the default relative tolerance is zero. Other KSP
types retain PETSc's default relative tolerance unless overridden through
its options database.

EDRIXS calls KSP's ``setFromOptions()`` after applying ``backend_kws``.
Corresponding PETSc options-database entries, such as ``ksp_atol``, ``ksp_rtol``,
``ksp_max_it``, and ``ksp_type``, take precedence over these settings.

``ksp_type`` accepts a PETSc KSP type such as ``'gmres'``; the available
algorithms and their tradeoffs are listed in the `KSP solver table
<https://petsc.org/release/manual/ksp/#tab-kspdefaults>`__.

To select GMRES instead of the default shifted PMINRES, pass::

    backend_kws={'ksp_type': 'gmres', 'linsys_tol': 1e-10}
