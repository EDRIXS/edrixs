===================
Quickstart tutorial
===================

If you are already familiar with multiplet calculations this page and our
:ref:`examples page <examples>` are a good part to start. Otherwise, checkout
our Pedagogical examples to learn the concepts
behind how edrixs works. The first example explains :ref:`Exact diagonalization 
<sphx_glr_auto_examples_example_0_ed_calculator.py>`.

Use edrixs as an ED calculator
------------------------------
edrixs can be used as a simple ED calculator to get eigenvalues (eigenvectors) of a many-body Hamiltonian with small size dimension (:math:`< 1,000`).
We will give an example to get eigenvalues for a :math:`t_{2g}`-orbital system (:math:`l_{eff}=1`). There are 6 orbitals including spin.

Launch your favorite python terminal::

    >>> import edrixs
    >>> import scipy
    >>> norb = 6
    >>> noccu = 2
    >>> Ud, JH = edrixs.UJ_to_UdJH(4, 1)
    >>> F0, F2, F4 = edrixs.UdJH_to_F0F2F4(Ud, JH)
    >>> umat = edrixs.get_umat_slater('t2g', F0, F2, F4)
    >>> emat = edrixs.atom_hsoc('t2g', 0.2)
    >>> basis = edrixs.get_fock_bin_by_N(norb, noccu)
    >>> H = edrixs.build_opers(4, umat, basis)
    >>> e1, v1 = scipy.linalg.eigh(H)
    >>> H += edrixs.build_opers(2, emat, basis)
    >>> e2, v2 = scipy.linalg.eigh(H)
    >>> print(e1)
    [1. 1. 1. 1. 1. 1. 1. 1. 1. 3. 3. 3. 3. 3. 6.]
    >>> print(e2)
    [0.890519 0.890519 0.890519 0.890519 0.890519 1.1      1.1      1.1
     1.183391 3.009481 3.009481 3.009481 3.009481 3.009481 6.016609]


Hello RIXS!
-----------

This is a "Hello World!" example for RIXS simulations at Ni (:math:`3d^8`) :math:`L_{2/3}` edges.
:math:`L_3` means transition from Ni-:math:`2p_{3/2}` to Ni-:math:`3d`, and
:math:`L_2` means transition from Ni-:math:`2p_{1/2}` to Ni-:math:`3d`.

.. plot:: pyplots/helloworld.py
   :include-source:
