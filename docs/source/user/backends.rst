.. _backends:

========
Backends
========

The backend specifies which linear algebra package is used in
EDRIXS. It changes how matrices are represented and which numerical algorithms
are used to solve a specific physical model.

The staged interface consists of :func:`edrixs.solvers.build_op` or
:func:`edrixs.solvers.get_ops`, followed by :func:`edrixs.solvers.ed`,
:func:`edrixs.solvers.xas`, and :func:`edrixs.solvers.rixs`.  Use ``backend``
to select an implementation and ``backend_kws`` to pass options specific to
that implementation and stage.

Available backends
------------------

``scipy``
   The fully implemented and recommended default.  Operators are SciPy CSR
   sparse matrices or ``LinearOperator`` objects.  It uses LOBPCG for ED,
   Lanczos continued fractions for XAS, and GMRES plus Lanczos for RIXS.

``dense``
   Constructs dense NumPy arrays and is intended as a compatibility option for
   small operators.  This backend is a work in progress: its ED, XAS, and RIXS
   entry points currently alias the SciPy iterative solvers rather than a
   separate exact dense eigensolver.

``petsc``
   Reserves an interface for PETSc matrices and vectors.  It is not yet
   implemented: calls require ``petsc4py`` and then raise
   ``NotImplementedError``.

The legacy functions whose names end in ``_fort`` use the compiled Fortran
workflows.  They predate the staged backend interface and are not selected
with ``backend="..."``.

Selecting a backend
-------------------

Operator construction defaults to ``backend="scipy"``.

For ``ed``, ``xas``, and ``rixs``, ``backend=None`` asks EDRIXS to infer the
backend from the supplied operators.  SciPy sparse matrices,
``LinearOperator`` objects, and NumPy arrays are accepted by the SciPy
backend.  PETSc objects are recognized when ``petsc4py`` is installed.  Mixed
or unrecognized operator types require an explicit, compatible backend.

SciPy backend options
---------------------

Operator construction
~~~~~~~~~~~~~~~~~~~~~

``tol`` (default ``1e-10``)
   Omit one- and two-body coefficients whose magnitude does not exceed this
   sparse-assembly threshold.  Accepted by ``build_op`` and ``get_ops``.

ED
~~

``num_evals`` (default ``1``)
   Number of lowest eigenpairs returned.  This is a direct argument to ``ed``,
   not an entry in ``backend_kws``.

``blocksize`` (default ``num_evals``)
   Number of eigenpairs requested internally from LOBPCG.  It must be at least
   ``num_evals``; extra results are discarded.

``tol`` (default ``1e-10``), ``maxiter`` (default ``200``)
   LOBPCG convergence tolerance and iteration limit.

``seed`` (default ``None``)
   Seed for the random initial block when ``initial_guess`` is absent.

``initial_guess`` (default ``None``)
   Starting block with shape ``(dimension, blocksize)``.

``suppress_lobpcg_warnings`` (default ``True``)
   Hide selected SciPy non-convergence warnings.  This does not improve
   convergence.

XAS
~~~

``nkryl`` (default ``200``)
   Maximum intermediate-state Lanczos dimension.

RIXS
~~~~

``nkryl`` (default ``200``)
   Maximum final-state Lanczos dimension.

``linsys_tol`` (default ``1e-9``)
   Relative GMRES stopping tolerance.  The absolute tolerance is zero on
   SciPy versions that distinguish it.

``linsys_maxiter`` (default ``50000``)
   Work limit for each GMRES solve.  Non-convergence raises an error.

``linsys_restart`` (default ``200``)
   GMRES Krylov-basis size before restart.  Larger values use more memory but
   can improve convergence.

.. _numerical-precision:

Numerical precision
-------------------

The SciPy backend relies on sparse iterative procedures.  Their output is
limited by stopping tolerances, Krylov dimensions, conditioning, and
floating-point arithmetic.  Successful completion does not imply that every
reported digit is exact.

Thermally populated states
~~~~~~~~~~~~~~~~~~~~~~~~~~

XAS and RIXS sum over every state supplied in ``eval_i`` and ``evec_i`` with
weights set up the Boltzman distribution at the specified temperature. The
User must be sure that enough states are included to represent the ground
state appropriately. In many cases, it's clear from the model that the ground
state must have a specific number of low-energy multiples. If not just compute
more ground state eigenvectors and ensure that the maximum energy one
substantially exceeds :math:`k_B T`.

LOBPCG eigenpairs
~~~~~~~~~~~~~~~~~

Set ``suppress_lobpcg_warnings=False`` while testing convergence.  Increase
``maxiter`` or ``blocksize`` when required, and compare runs with different
initial blocks.  Check every eigenpair using

.. math::

   r_i = \lVert H v_i-E_i v_i\rVert,

and verify that ``evec_i.conj().T @ evec_i`` is close to the identity.  A
larger block is often useful around degeneracies, but only ``num_evals``
states are returned and included in later thermal averages.

Lanczos and re-orthogonalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In exact arithmetic Lanczos vectors are orthogonal.  Roundoff can gradually
destroy that orthogonality, producing repeated Ritz values or inaccurate
spectral weight in long recurrences or near degeneracies.  The current EDRIXS
Lanczos implementation does not perform explicit re-orthogonalization and
does not expose a re-orthogonalization option.

Increase ``nkryl`` in steps until peak positions and weights are stable.
Increasing it indefinitely is not guaranteed to help because truncation error
falls while loss of orthogonality can grow.  Sharper spectra and smaller
broadening generally require a larger Krylov space.

GMRES in RIXS
~~~~~~~~~~~~~

RIXS solves a correction-vector system for every incident energy and retained
state.  Tighten ``linsys_tol`` until the spectrum is stable.  If GMRES fails,
increase ``linsys_maxiter`` and consider a larger ``linsys_restart``.  A
smaller ``gamma_c`` can make the system harder to solve, but ``gamma_c`` is a
physical lifetime broadening and should not be changed merely to hide
non-convergence.

Convergence checklist
~~~~~~~~~~~~~~~~~~~~~

* Keep all thermally relevant states and complete degenerate manifolds.
* Inspect ED residuals and SciPy convergence warnings.
* Converge ``nkryl`` for both XAS and RIXS.
* Converge the RIXS GMRES settings independently.
* Check the operator-assembly ``tol`` when small matrix elements matter.
* Use energy grids fine enough to resolve features set by ``gamma_c`` and
  ``gamma_f``; these broadenings are physical parameters, not tolerances.
* Report solver settings and changes in observables, rather than relying on
  the number of digits written to an output file.
