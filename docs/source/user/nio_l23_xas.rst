.. _nio-l23-xas:

**********************************************
NiO crystal-field :math:`L_{2,3}` XAS & RIXS
**********************************************

The calculation treats NiO as
an ionic Ni :math:`3d^8` site in a cubic field and excites a :math:`2p`
electron into the :math:`3d`. It includes an effective exchange
field acting along the :math:`<112>` direction.

The RIXS spectrum is calculated at the global maximum of the XAS and summed
over two orthogonal outgoing linear polarizations.

For RIXS we compute the energy loss spectrum for an incident energy
matching the maximum of the XAS.

.. plot:: pyplots/nio_l23_xas.py
   :context: reset
   :include-source: True
