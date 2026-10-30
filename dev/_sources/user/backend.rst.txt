.. _backend:

*******
Backend
*******

The :ref:`introduction` describes the four stages of a calculation.  Stages two
to four -- building operators, diagonalizing, and assembling spectra --
delegate their numerical work to a backend selected with ``backend`` when
calling ``build_op`` or ``get_ops``:

* ``'dense'`` builds full NumPy matrices and uses exact dense
  diagonalization and explicit eigenstate sums.  It is intended for small
  problems, reference calculations, and pedagogy.
* ``'scipy'`` builds SciPy CSR matrices and uses LOBPCG, Lanczos, and GMRES.
  It is the default backend.
* ``'fortran'`` uses the native MPI solvers and communicates through files in
  the current working directory.  We provide the function
  :func:`edrixs.get_ops_disk <edrixs.solvers.get_ops_disk>` to load the setup
  from disk if desired.
* ``'petsc'`` uses distributed PETSc matrices, SLEPc eigensolvers, and PETSc
  Krylov solvers.  It is intended for larger Hilbert spaces and MPI runs.

Backend-specific numerical controls belong in a ``backend_kws`` mapping. The
mapping is specific to both the backend *and the operation*.

.. _lazy-operators:

Operator factories and solver dispatch
======================================

``build_op`` returns a callable that constructs the matrix when called. This
lets callers defer matrix allocation and recreate matrices when needed. Its
``.backend`` attribute identifies the numerical backend and the matrix
representation it will produce.

Numerical ``get_ops`` returns these callables for both Hamiltonians
and every transition component. Public ``ed`` and ``rixs`` dispatch using
``hmat_i.backend``; ``xas`` uses ``hmat_n.backend``. They pass the factories
unchanged to the selected backend and assume all operators are suitable for it.
Dense, SciPy, and PETSc ED call every supplied factory at calculation
start. Repeated arguments may be instantiated independently. PETSc XAS and RIXS
accept ``matrix_lifetime`` in ``backend_kws``; see :doc:`backends/petsc`.
All staged solver entry points, including backend-specific functions, require
factories or callable Fortran disk handles.

For matrix arithmetic or direct access to matrix methods, call the factory
explicitly::

    number_op = edrixs.build_op(number_coefficients, None, basis, backend='dense')
    number_matrix = number_op()

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
* :doc:`PETSc keywords <backends/petsc>`

These references apply to the staged ``get_ops`` / ``ed`` / ``xas`` / ``rixs``
interface.  The legacy ``ed_*_fort``, ``xas_*_fort``, and ``rixs_*_fort``
functions retain their existing explicit parameters.

.. toctree::
   :hidden:

   backends/scipy
   backends/dense
   backends/fortran
   backends/petsc
