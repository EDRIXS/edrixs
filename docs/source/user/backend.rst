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

Preparing and using operators
============================

Choose the backend when preparing operators, then pass the returned objects
directly to the solvers::

    hmat_i, hmat_n, trans_ops = edrixs.get_ops(*problem, backend='scipy')
    eval_i, evec_i = edrixs.ed(hmat_i, shift=shift, num_evals=3)

Here ``problem`` and ``shift`` come from a model constructor as shown in
:ref:`physical-model`. You do not repeat the backend choice in ``ed``, ``xas``,
or ``rixs``. The operators carry that choice in their ``.backend`` labels;
prepare all operators for a calculation with the same backend.

For dense, SciPy, and PETSc calculations, these objects are instructions for
building a matrix, rather than matrices already occupying memory. In Python
they are callables (functions you can call with no arguments), also called
*operator factories*. The solvers call them when needed, and each call builds
a fresh matrix. This lets the backend release large matrices and reconstruct
them later. PETSc offers a memory-versus-computation choice through
``matrix_lifetime``; see :doc:`backends/petsc`.

Fortran instead uses handles to solver input files. Pass those handles to the
same solver functions; see :doc:`backends/fortran` for its file-based workflow.

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

Backend settings reference
==========================

Choose the backend whose solver and matrix settings you need:

* :doc:`backends/scipy`
* :doc:`backends/dense`
* :doc:`backends/fortran`
* :doc:`backends/petsc`

These references apply to the staged ``get_ops`` / ``ed`` / ``xas`` / ``rixs``
interface.  The legacy ``ed_*_fort``, ``xas_*_fort``, and ``rixs_*_fort``
functions retain their existing explicit parameters.

.. toctree::
   :hidden:

   backends/scipy
   backends/dense
   backends/fortran
   backends/petsc
