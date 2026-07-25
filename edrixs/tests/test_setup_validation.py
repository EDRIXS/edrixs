"""Unit tests for model-definition validation in the solver-neutral setup layer.

The ``setup_*`` functions define orbital spaces, interactions, and photon
transition targets before either dense or SciPy sparse operators are built.
These tests check invalid model choices at that earliest workflow stage.
"""

import pytest

from edrixs.solvers import setup_2v1c, setup_siam


def test_setup_2v1c_rejects_unknown_transition_target():
    """Reject an undefined core-to-valence target in ``setup_2v1c``.

    The selected target determines where the photon transition operator is
    placed before ``ops`` builds dense or sparse many-electron operators.  This
    test stops an ambiguous model before any ED, XAS, or RIXS calculation.
    """
    with pytest.raises(Exception, match="Unknown trans_to_which"):
        setup_2v1c(("s", "s", "p"), trans_to_which=3)


def test_setup_siam_rejects_unknown_type():
    """Reject an unsupported SIAM representation during model setup.

    ``setup_siam`` must know whether impurity, bath, and hybridization data are
    supplied in the supported compact or full form before ``ops`` can construct
    Hamiltonians.  This test catches an invalid choice at that boundary.
    """
    with pytest.raises(Exception, match="Unknown siam_type"):
        setup_siam(("s", "p"), 1, siam_type=7)
