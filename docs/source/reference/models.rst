models
======

Model constructors define the physical problem without building many-body
operators or running solvers. They return eight values: seven inputs describing
the Hamiltonians, Fock bases, and transition operators, plus a scalar ``shift``
equal to the filled core's shell-level energy.

Pass the first seven outputs to ``get_ops``. For dense, SciPy, or PETSc
diagonalization, optionally pass ``shift`` to ``ed`` to reduce numerical errors.
The returned eigenvalues retain their original energy reference, so the shift
does not change the physics. Fortran ED ignores a nonzero shift with a warning.

.. automodule:: edrixs.models
   :members:
