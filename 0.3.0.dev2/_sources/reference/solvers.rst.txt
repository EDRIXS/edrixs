solvers
=======

``build_op`` returns a zero-argument ``functools.partial`` with a ``.backend``
attribute. Numerical ``get_ops`` returns these factories for both Hamiltonians
and each transition component. ``ed`` and ``rixs`` dispatch using
``hmat_i.backend``; ``xas`` uses ``hmat_n.backend``. The backend receives the
factories unchanged. All staged backend solvers require factories.
Dense, SciPy, and PETSc ED construct all inputs at calculation start.
PETSc XAS and RIXS accept ``matrix_lifetime`` in ``backend_kws``.
Fortran disk handles are callable, returning themselves; native Fortran manages
actual matrix construction.

Call a numerical factory explicitly to obtain a concrete matrix for arithmetic.

.. automodule:: edrixs.solvers
   :members:
