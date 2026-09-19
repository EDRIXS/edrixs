.. _nio-siam-xas-rixs:

*******************************************
NiO Anderson impurity model XAS and RIXS
*******************************************

This example uses :func:`~edrixs.models.model_siam` to calculate Ni
:math:`L_{2,3}`-edge XAS and RIXS for a NiO Anderson impurity model. The model
contains a correlated Ni :math:`3d` shell and one bath site representing ten
symmetry-adapted O :math:`2p` ligand spin-orbitals. Hybridization allows the
calculation to include ligand-to-metal charge-transfer configurations.

.. plot:: pyplots/nio_siam_xas_rixs.py
   :context: reset
   :include-source: True
