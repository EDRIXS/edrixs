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

Options are specific to the stage where they are supplied.

``build_op`` and ``get_ops``
============================

.. backend-options:: petsc get_ops

``nnz_guess_per_row`` controls PETSc matrix preallocation.  Leaving it unset
uses an estimate from the retained one- and two-body coefficients.  A
``mat_type`` such as ``'aijcusparse'`` can select an accelerated matrix format;
see PETSc's `matrix type overview
<https://petsc.org/release/manual/mat/#basic-matrix-operations>`__.

``ed``
======

Pass the final model output as ``ed(hmat_i, shift=shift)`` to subtract the
filled-core shell energy from a copy of the distributed Hamiltonian before
diagonalization. Returned eigenvalues retain the original energy reference.
With a nonzero ``shift``, ``eigval_tol`` is an absolute residual tolerance for
the shifted problem, avoiding division by a near-zero centered eigenvalue.
With ``shift=0``, it is a relative residual tolerance, scaled by the eigenvalue
magnitude.

.. backend-options:: petsc ed

The eigensolver uses SLEPc ``EPS`` with a Hermitian problem type and requests
the smallest-real eigenpairs.  ``ncv`` sets the Krylov subspace size described
by `EPS.setDimensions
<https://slepc.upv.es/release/slepc4py/reference/slepc4py.SLEPc.EPS.html#slepc4py.SLEPc.EPS.setDimensions>`__.

``xas``
=======

XAS and RIXS center distributed matrix copies using the lowest retained initial
energy and the first incident energy, which should lie near resonance. RIXS
also centers the final Hamiltonian. Spectra are evaluated on these relative
energies before restoring the original reference in returned pole dictionaries.
Input matrices are unchanged.

.. backend-options:: petsc xas

``rixs``
========

.. backend-options:: petsc rixs

``linsys_tol`` sets KSP's absolute tolerance, ``atol``. The relative tolerance
remains at PETSc's default unless overridden through its options database;
the solver can therefore converge through the relative criterion before
reaching the absolute tolerance.

EDRIXS calls KSP's ``setFromOptions()`` after applying ``backend_kws``.
Corresponding PETSc options-database entries, such as ``ksp_atol``, ``ksp_rtol``,
``ksp_max_it``, and ``ksp_type``, take precedence over these settings.

``ksp_type`` accepts a PETSc KSP type such as ``'gmres'``; the available
algorithms and their tradeoffs are listed in the `KSP solver table
<https://petsc.org/release/manual/ksp/#tab-kspdefaults>`__.
