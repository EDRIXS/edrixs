.. _backend:

*******
Backend
*******

The :ref:`introduction` describes the four stages of a calculation.  Stages two
to four -- building operators, diagonalizing, and assembling spectra -- do not
implement their own linear algebra.  They delegate it to a backend selected
with the ``backend`` keyword, so one argument takes the same script from a toy
problem to a production one:

* ``'dense'`` -- full NumPy matrices.  Use it for small problems
  (dimension :math:`\lesssim 1000`) and for pedagogy.
* ``'scipy'`` -- SciPy sparse matrices with Lanczos / Krylov solvers.  The
  default, and a good choice for single-site and small cluster models.
* ``'petsc'`` -- distributed PETSc/SLEPc sparse matrices.  Needed for the large
  Hilbert spaces of Anderson impurity models; runs in parallel when the script
  is launched with ``mpirun``.

Purely numerical knobs (Krylov dimension, solver tolerances, ...) are passed per
backend through a ``backend_kws`` dictionary, so that the physical arguments to
``ed``, ``xas`` and ``rixs`` stay backend independent.
