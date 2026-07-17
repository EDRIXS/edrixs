#!/usr/bin/env python

import numpy as np
import collections
import edrixs
from mpi4py import MPI

if __name__ == "__main__":
    '''
    U-5f2 :math:`L_{3}`-edge, :math:`2p_{3/2}\\rightarrow 6d` transition.

    How to run
    ----------
    mpiexec -n 4 python run_rixs_fsolver.py
    '''
    # Occupancy of U 5f orbitals
    noccu = 2
    res = edrixs.get_atom_data('U', v_name=('5f', '6d'), v_noccu=(noccu, 0),
                               edge='L3', trans_to_which=2, label=('f', 'd', 'p'))

    slat_i = collections.OrderedDict(res['slater_i'])
    slat_n = collections.OrderedDict(res['slater_n'])

    # Initial Hamiltonian
    slat_i['F2_ff'] = slat_i['F2_ff'] * 0.77
    slat_i['F4_ff'] = slat_i['F4_ff'] * 0.77
    slat_i['F6_ff'] = slat_i['F6_ff'] * 0.77
    slat_i['F0_ff'] = edrixs.get_F0('f', slat_i['F2_ff'], slat_i['F4_ff'], slat_i['F6_ff'])

    # Intermediate Hamiltonian
    slat_n['F2_ff'] = slat_n['F2_ff'] * 0.77
    slat_n['F4_ff'] = slat_n['F4_ff'] * 0.77
    slat_n['F6_ff'] = slat_n['F6_ff'] * 0.77

    slat_n['F0_ff'] = edrixs.get_F0('f', slat_n['F2_ff'], slat_n['F4_ff'], slat_n['F6_ff'])
    slat_n['F0_fd'] = edrixs.get_F0('fd', slat_n['G1_fd'], slat_n['G3_fd'], slat_n['G5_fd'])
    slat_n['F0_fp'] = edrixs.get_F0('fp', slat_n['G2_fp'], slat_n['G4_fp'])
    slat_n['F0_dp'] = edrixs.get_F0('dp', slat_n['G1_dp'], slat_n['G3_dp'])

    slater = (list(slat_i.values()), list(slat_n.values()))

    # Spin-Orbit Coupling (SOC) zeta
    zeta_f_i = res['v_soc_i'][0]
    zeta_f_n = res['v_soc_n'][0]

    zeta_d_i = res['v_soc_i'][1]
    zeta_d_n = res['v_soc_n'][1]

    # RIXS settings
    thin, thout, phi = 45 / 180.0 * np.pi, 45 / 180.0 * np.pi, 0.0
    gamma_c = res['gamma_c'][0]
    gamma_f = 0.1
    ominc_xas = np.linspace(-10, 20, 1000)
    ominc_rixs = np.linspace(0, 10, 10)
    eloss = np.linspace(-0.2, 5, 1000)

    poltype_xas = [('isotropic', 0.0)]

    poltype_rixs = [('linear', 0.0, 'linear', 0.0),
                    ('linear', 0.0, 'linear', np.pi / 2.0)]

    shell_name = ('f', 'd', 'p32')

    # mpi4py env
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Run ED
    eval_i, denmat = edrixs.ed_2v1c_fort(
        comm, shell_name, shell_level=(0, 5.0, 0),
        v1_soc=(zeta_f_i, zeta_f_n), v2_soc=(zeta_d_i, zeta_d_n),
        v_tot_noccu=noccu, slater=slater, ed_solver=2, neval=20,
        nvector=9, ncv=50, idump=True
    )
    if rank == 0:
        print('eigvals:', eval_i)
        print('occupancy numbers:', denmat[0].diagonal())

    # Run XAS
    xas, xas_poles = edrixs.xas_2v1c_fort(
        comm, shell_name, ominc_xas, gamma_c=gamma_c,
        v_tot_noccu=noccu, trans_to_which=2, thin=thin, phi=phi,
        pol_type=poltype_xas, num_gs=9, nkryl=200, temperature=300
    )
    if rank == 0:
        np.savetxt('xas.dat', np.concatenate((np.array([ominc_xas]).T, xas), axis=1))
        edrixs.dump_poles(xas_poles, 'xas_poles')

    # Run RIXS
    rixs, rixs_poles = edrixs.rixs_2v1c_fort(
        comm, shell_name, ominc_rixs, eloss, gamma_c=gamma_c, gamma_f=gamma_f,
        v_tot_noccu=noccu, trans_to_which=2, thin=thin, thout=thout, phi=phi,
        pol_type=poltype_rixs, num_gs=9, nkryl=200, linsys_max=1000, temperature=300
    )

    if rank == 0:
        edrixs.dump_poles(rixs_poles, 'rixs_poles')
        rixs_pi = np.sum(rixs[:, :, 0:2], axis=2)
        np.savetxt('rixs_pi.dat', np.concatenate((np.array([eloss]).T, rixs_pi.T), axis=1))
        edrixs.plot_rixs_map(rixs_pi, ominc_rixs, eloss, "rixsmap_pi.pdf")
