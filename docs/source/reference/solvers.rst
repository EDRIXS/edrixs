solvers
=======

``build_op`` returns a zero-argument ``functools.partial`` with a ``.backend``
attribute. Numerical ``get_ops`` returns these factories for both Hamiltonians
and each transition component. ``ed``, ``xas``, and ``rixs`` take the factories,
validate their labels, and construct their matrices. Fortran disk handles carry the same label and are passed
through without being called.

Call a numerical factory explicitly to obtain a concrete matrix for arithmetic.
See :ref:`lazy-operators` for examples, lifetime semantics, and migration notes.

.. automodule:: edrixs.solvers
   :members:
