models
======
Model constructors return eight values: seven orbital-space problem inputs
and a scalar ``shift`` equal to the filled core's shell-level energy. Coulomb
contributions are excluded from this shift, and the integrals remain unshifted.
Use it explicitly for dense, SciPy, or PETSc diagonalization::

    *problem, shift = edrixs.model_1v1c(('p', 's'), shell_level=(0, -100))
    hmat_i, hmat_n, trans_ops = edrixs.get_ops(*problem, backend='scipy')
    eval_i, evec_i = edrixs.ed(hmat_i, shift=shift)

The returned eigenvalues retain their absolute energy reference. Fortran
callers pass only ``problem`` to ``get_ops`` and omit ``shift`` from ``ed``.

.. automodule:: edrixs.models
   :members:
