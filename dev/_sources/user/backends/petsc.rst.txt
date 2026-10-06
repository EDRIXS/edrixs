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

The EDRIXS PETSc backend is currently an unimplemented contract.  It recognizes
PETSc ``Mat`` and ``Vec`` objects, but its operations raise
``NotImplementedError`` and accept no backend-specific options.  The keyword
contract will be documented here when those operations are implemented.
