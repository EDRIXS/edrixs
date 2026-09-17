__all__ = ['write_spectrum', 'plot_rixs_map']

import numpy as np
import matplotlib.pyplot as plt

from . import poles


def write_spectrum(file_list, omega_mesh, gamma_mesh, T=1.0, fname='spectrum.dat',
                   om_shift=0.0, fmt_float='{:.15f}'):
    """
    Reading poles :math:`\\alpha` and :math:`\\beta`, and calculate
    the spectrum using continued fraction formula,

    .. math::
        I(\\omega_{i}) =-\\frac{1}{\\pi}\\text{Im} \\left[ \\frac{1}{x - \\alpha_{0} -
        \\frac{\\beta_{1}^2}{x-\\alpha_{1} - \\frac{\\beta_{2}^2}{x-\\alpha_{2} - ...}} }\\right],

    where, :math:`x = \\omega_{i}+i\\Gamma_{i} + E_{g}`.

    Parameters
    ----------
    file_list: list of string
        Name of poles file.
    omega_mesh: 1d float array
        The frequency mesh.
    gamma_mesh: 1d float array
        The broadening factor, in general, it is frequency dependent.
    T: float (default: 1.0K)
        Temperature (K).
    fname: str (default: 'spectrum.dat')
        File name to store spectrum.
    om_shift: float (default: 0.0)
        Energy shift.
    fmt_float: str (default: '{:.15f}')
        Format for printing float numbers.
    """

    from .fortran_backend.isostream_fortran import read_poles_from_file

    pole_dict = read_poles_from_file(file_list)
    spectrum = poles.get_spectra_from_poles(pole_dict, omega_mesh, gamma_mesh, T)

    space = "  "
    fmt_string = (fmt_float + space) * 2 + '\n'
    f = open(fname, 'w')
    for i in range(len(omega_mesh)):
        f.write(fmt_string.format(omega_mesh[i] + om_shift, spectrum[i]))
    f.close()


def plot_rixs_map(rixs_data, ominc_mesh, eloss_mesh, fname='rixsmap.pdf'):
    """
    Given 2d RIXS data, plot a RIXS map and save it to a pdf file.

    Parameters
    ----------
    rixs_data: 2d float array
        Calculated RIXS data as a function of incident energy and energy loss.
    ominc_mesh: 1d float array
        Incident energy mesh.
    eloss_mesh: 1d float array
        Energy loss mesh.
    fname: string
        File name to save RIXS map.
    """

    fig, ax = plt.subplots()
    a, b, c, d = min(eloss_mesh), max(eloss_mesh), min(ominc_mesh), max(ominc_mesh)
    m, n = np.array(rixs_data).shape
    if len(ominc_mesh) == m and len(eloss_mesh) == n:
        plt.imshow(
            rixs_data, extent=[a, b, c, d], origin='lower', aspect='auto',
            cmap='rainbow', interpolation='gaussian'
        )
        plt.xlabel(r'Energy loss (eV)')
        plt.ylabel(r'Energy of incident photon (eV)')
    elif len(eloss_mesh) == m and len(ominc_mesh) == n:
        plt.imshow(
            rixs_data, extent=[c, d, a, b], origin='lower', aspect='auto',
            cmap='rainbow', interpolation='gaussian'
        )
        plt.ylabel(r'Energy loss (eV)')
        plt.xlabel(r'Energy of incident photon (eV)')
    else:
        raise Exception(
            "Dimension of rixs_data is not consistent with ominc_mesh or eloss_mesh"
        )

    plt.savefig(fname)
