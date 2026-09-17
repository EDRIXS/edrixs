__all__ = ['get_spectra_from_poles', 'merge_pole_dicts']

import numpy as np

from .utils import boltz_dist


def get_spectra_from_poles(poles_dict, omega_mesh, gamma_mesh, temperature):
    """
    Given the dict of poles, calculate XAS or RIXS spectra using continued fraction formula,

    .. math::
        I(\\omega_{i}) =-\\frac{1}{\\pi}\\text{Im} \\left[ \\frac{1}{x - \\alpha_{0} -
        \\frac{\\beta_{1}^2}{x-\\alpha_{1} - \\frac{\\beta_{2}^2}{x-\\alpha_{2} - ...}} }\\right],

    where, :math:`x = \\omega_{i}+i\\Gamma_{i} + E_{g}`.

    Parameters
    ----------
    poles_dict: dict
        Dict containing information of poles, which are calculated from
        xas_fsolver and rixs_fsolver.
        This dict is constructed by
        :func:`fortran_backend.isostream_fortran.read_poles_from_file`.
    omega_mesh: 1d float array
        Energy grid.
    gamma_mesh: 1d float array
        Life-time broadening.
    temperature: float number
        Temperature (K) for boltzmann distribution.

    Returns
    -------
    spectra: 1d float array
        The calculated XAS or RIXS spectra.

    See also
    --------
    fortran_backend.isostream_fortran.read_poles_from_file:
        Read XAS or RIXS poles files.
    """
    nom = len(omega_mesh)
    spectra = np.zeros(nom, dtype=np.float64)
    gs_dist = boltz_dist(poles_dict['eigval'], temperature)
    ngs = len(poles_dict['eigval'])
    for i in range(ngs):
        tmp_vec = np.zeros(nom, dtype=complex)
        neff = poles_dict['npoles'][i]
        alpha = poles_dict['alpha'][i]
        beta = poles_dict['beta'][i]
        eigval = poles_dict['eigval'][i]
        norm = poles_dict['norm'][i]
        for j in range(neff-1, 0, -1):
            tmp_vec = (
                beta[j-1]**2 / (omega_mesh + 1j * gamma_mesh + eigval - alpha[j] - tmp_vec)
            )
        tmp_vec = (
            1.0 / (omega_mesh + 1j * gamma_mesh + eigval - alpha[0] - tmp_vec)
        )
        spectra[:] += -1.0 / np.pi * np.imag(tmp_vec) * norm * gs_dist[i]

    return spectra


def merge_pole_dicts(list_pole_dict):
    """
    Given a list of dict of poles, merge them into one dict of poles

    Parameters
    ----------
    list_pole_dict:  list of dict
        Dict containing information of poles, which are calculated from
        xas_fsolver and rixs_fsolver.

    Returns
    -------
    new_pole_dict: dict of poles
        New dict of poles.
    """
    new_pole_dict = {
        'eigval': [],
        'npoles': [],
        'norm': [],
        'alpha': [],
        'beta': []
    }
    for poles_dict in list(list_pole_dict):
        new_pole_dict['eigval'].extend(poles_dict['eigval'])
        new_pole_dict['npoles'].extend(poles_dict['npoles'])
        new_pole_dict['norm'].extend(poles_dict['norm'])
        new_pole_dict['alpha'].extend(poles_dict['alpha'])
        new_pole_dict['beta'].extend(poles_dict['beta'])

    return new_pole_dict
