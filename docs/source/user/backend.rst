.. _backend:

*******
Backend
*******

The :ref:`introduction` describes the four stages of a calculation.  Stages two
to four -- building operators, diagonalizing, and assembling spectra --
delegate their numerical work to a backend selected with ``backend``:

* ``'dense'`` returns full NumPy matrices and uses exact dense
  diagonalization and explicit eigenstate sums.  It is intended for small
  problems, reference calculations, and pedagogy.
* ``'scipy'`` returns SciPy sparse matrices and uses LOBPCG, Lanczos, and GMRES.
  It is the default backend.
* ``'fortran'`` uses the native MPI solvers and communicates through files in
  the current working directory.  We provide the function
  :func:`edrixs.get_ops_disk <edrixs.solvers.get_ops_disk>` to load the setup
  from disk if desired.
* ``'petsc'`` reserves the PETSc backend contract but is not yet implemented.

Backend-specific numerical controls belong in a ``backend_kws`` mapping.  The
mapping is specific to both the backend *and the operation*.

Numerical precision
===================

XAS and RIXS thermally average over the retained initial eigenstates.  Request
enough eigenpairs from ``ed`` to keep every thermally populated state in the
ground-state manifold, including complete degenerate or nearly degenerate
multiplets; otherwise the Boltzmann average is incomplete.

EDRIXS uses established numerical linear-algebra packages, but floating-point
and iterative solvers are not exact.  If numerical accuracy is a concern,
repeat the calculation with varied solver settings and check that the result
is unchanged.  Relevant settings include ``tol`` or ``eigval_tol``,
``maxiter``, ``blocksize`` or ``ncv``, ``nkryl``, and the RIXS ``linsys_*``
options documented below.

.. _backend-options:

Backend keyword reference
=========================

Choose the backend whose keywords you need:

* :doc:`SciPy keywords <backends/scipy>`
* :doc:`Dense keywords <backends/dense>`
* :doc:`Fortran keywords <backends/fortran>`
* :doc:`PETSc status and keywords <backends/petsc>`

These references apply to the staged ``get_ops`` / ``ed`` / ``xas`` / ``rixs``
interface.  The legacy ``ed_*_fort``, ``xas_*_fort``, and ``rixs_*_fort``
functions retain their existing explicit parameters.

.. toctree::
   :hidden:

   backends/scipy
   backends/dense
   backends/fortran
   backends/petsc
