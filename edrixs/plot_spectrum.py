__all__ = ['plot_rixs_map']

import numpy as np
import matplotlib.pyplot as plt


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
