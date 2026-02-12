==================
Expectation values
==================

Here we cover computing expectation values. For small problems,
this can be done using pure python, as explained in the Pedagogical examples.
EDRIXS also provides methods to do this with the Fortran solver, which can handle
problems with more orbitals. Here, we choose to repeat the simple and computationally
inexpensive example of computing the :math:`d`-electron occupation in NiO. This is
done using the ``opavg.x`` executable provided by EDRIXS, which reads files
from the working directory and writes a file descrbing the desired eigenvalues.

Create eigenvectors
-------------------
The first step is to run :code:`edrixs.ed_siam_fort`, as explained in
the :ref:`Pedagogical example for NiO
<sphx_glr_auto_examples_example_3_AIM_XAS.py>`. This will write files describing the
eigenvectors as ``eigenvec.1``, ``eigenvec.2``, ``eigenvec.3``, ``...``. Copy these
files to a separate folder called, for example, ``opavg``. Change the working directory
to this location. The rest of the example will occur in this directory.
If you would like to obtain all the eigenvectors for NiO, pass the keyword argument
:code:`nvector=190` to :code:`edrixs.ed_siam_fort`. 

Define Fock basis
-----------------
:code:`opavg.x` needs to know the Fock basis for our NiO example, which has 20 spin
orbitals hosting 18 electron (8 Ni electrons and 10 electrons in the ligand bath). This
is done by creating the file ``fock.in`` file as follows.

.. code-block:: python

   import edrixs
   v_noccu = 18
   ntot = 20
   edrixs.write_fock_dec_by_N(ntot, v_noccu, 'fock_i.in')

Define operators
----------------
The operators computed by :code:`opavg.x` are defined by ``hopping_i.in`` and
``coulomb_i.in``. The number operator involves only two fermion terms in
``hopping_i.in``, so we can set ``coulomb_i.in`` to zero, taking care to
define an appropriately sized matrix. In the real harmonic basis, we can
count the occupation of orbitals using a identity matrix, which we then
transform to the default complex harmonic basis.

.. code-block:: python

   import edrixs
   import numpy as np
   ndorb = 10
   ntot = 20
   nd_real_harmoic_basis = np.identity(ndorb, dtype=complex)
   nd_complex_harmoic_basis = edrixs.cb_op(nd_real_harmoic_basis,
                                           edrixs.tmat_r2c('d', True))
   emat = np.zeros((ntot, ntot), dtype=complex)
   emat[:ndorb, :ndorb] = nd_complex_harmoic_basis
   edrixs.write_emat(emat, 'hopping_i.in')
   
   umat = np.zeros((ntot, ntot, ntot, ntot), dtype=complex)
   edrixs.write_umat(umat, 'coulomb_i.in')


Setup configuration
-------------------
The diagonalization routine is controled by creating a file ``config.in`` with
the following contents:

.. code-block:: console

   &control
   ed_solver     =   2
   num_val_orbs  =   20
   ncv           =   50
   neval         =   190
   nvector       =   190
   num_gs        =   190
   maxiter       =   1000
   eigval_tol    =   1E-10
   idump         =   .true.
   &end

See :code:`edrixs.ed_siam_fort` for an explanation of the meaning of
these parameters.

Run and examine values
----------------------
Executing :code:`opavg.x` from the terminal will generate
:code:`opavg.dat` containing the required eigenvalues.