import numpy as np

from edrixs.wannier_ham import HR


def test_hr_from_file(tmp_path):
    wannier_file = tmp_path / "wannier90_hr.dat"
    wannier_file.write_text(
        "Generated for testing\n"
        "1\n"
        "1\n"
        "1\n"
        "0 0 0 1 1 2.0 -0.5\n"
    )

    result = HR.from_file(wannier_file)

    assert result.nwann == 1
    assert result.nrpt == 1
    np.testing.assert_array_equal(result.deg_rpt, [1])
    np.testing.assert_allclose(result.hr[0, 0, 0], 2.0 - 0.5j)
