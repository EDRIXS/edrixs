"""Calculate NiO Anderson-impurity-model XAS and RIXS spectra."""

import matplotlib.pyplot as plt
import numpy as np

import edrixs


# Model size: a nominal d8 impurity plus one filled ten-orbital bath.
nd = 8
norb_d = 10
nbath = 1
v_noccu = nd + nbath * norb_d
shell_name = ('d', 'p')

# Screened Slater integrals and monopole interactions from example 03.
info = edrixs.get_atom_data('Ni', '3d', nd, edge='L3')
F2_dd = 0.8 * info['slater_i'][1][1]
F4_dd = 0.8 * info['slater_i'][2][1]
U_dd = 7.3
F0_dd = U_dd + edrixs.get_F0('d', F2_dd, F4_dd)

F2_dp = 0.8 * info['slater_n'][4][1]
G1_dp = 0.8 * info['slater_n'][5][1]
G3_dp = 0.8 * info['slater_n'][6][1]
U_dp = 8.5
F0_dp = U_dp + edrixs.get_F0('dp', G1_dp, G3_dp)
slater = (
    [F0_dd, F2_dd, F4_dd],
    [F0_dd, F2_dd, F4_dd, F0_dp, F2_dp, G1_dp, G3_dp],
)

# Convert the charge-transfer energy into impurity, ligand, and core levels.
Delta = 4.7
E_d, E_L = edrixs.CT_imp_bath(U_dd, Delta, nd)
E_dc, E_Lc, E_p = edrixs.CT_imp_bath_core_hole(
    U_dd, U_dp, Delta, nd
)

# Work in the real-harmonic basis used for the bath parameterization.
trans_c2n = edrixs.tmat_c2r('d', True)
ten_dq = 0.56
orbital_energies = np.repeat(
    [0.6 * ten_dq, -0.4 * ten_dq, -0.4 * ten_dq,
     0.6 * ten_dq, -0.4 * ten_dq],
    2,
)
crystal_field = np.diag(orbital_energies).astype(complex)
soc = edrixs.cb_op(
    edrixs.atom_hsoc('d', info['v_soc_i'][0]), trans_c2n
)
imp_mat = crystal_field + soc + E_d * np.eye(norb_d)
imp_mat_n = crystal_field + soc + E_dc * np.eye(norb_d)

# Bath crystal field and impurity--bath hybridization.
ten_dq_bath = 1.44
bath_splitting = np.repeat(
    [0.6 * ten_dq_bath, -0.4 * ten_dq_bath, -0.4 * ten_dq_bath,
     0.6 * ten_dq_bath, -0.4 * ten_dq_bath],
    2,
)
bath_level = E_L + bath_splitting[np.newaxis, :]
bath_level_n = E_Lc + bath_splitting[np.newaxis, :]

hyb = np.zeros((nbath, norb_d), dtype=complex)
hyb[0] = np.repeat([2.06, 1.21, 1.21, 2.06, 1.21], 2)

# The core-level shift places the calculated spectrum near the Ni L edge.
edge_shift = 857.6
c_level = -edge_shift - 5 * E_p
problem = edrixs.model_siam(
    shell_name,
    nbath,
    siam_type=0,
    v_noccu=v_noccu,
    c_level=c_level,
    c_soc=info['c_soc'],
    trans_c2n=trans_c2n,
    imp_mat=imp_mat,
    imp_mat_n=imp_mat_n,
    bath_level=bath_level,
    bath_level_n=bath_level_n,
    hyb=hyb,
    slater=slater,
    ext_B=np.array([0.0, 0.0, 0.12]),
    on_which='spin',
    sparse_U=True,
)

backend = 'scipy'
hmat_i, hmat_n, trans_ops = edrixs.get_ops(*problem, backend=backend)
eval_i, evec_i = edrixs.ed(
    hmat_i,
    num_evals=3,
    backend=backend,
    backend_kws={'tol': 1e-10, 'maxiter': 1000, 'seed': 0},
)

# Isotropic XAS with the geometry, temperature, and broadening of example 03.
temperature = 300
thin = 0.0
phi = 0.0
energy = edge_shift + np.linspace(-15, 25, 1000)
gamma_c = 0.48 / 2
xas = edrixs.xas(
    eval_i, evec_i, hmat_n, trans_ops, energy,
    gamma_c=gamma_c,
    thin=thin,
    phi=phi,
    pol_type=[('isotropic', 0)],
    temperature=temperature,
    backend=backend,
    backend_kws={'nkryl': 120},
)[:, 0]
xas /= xas.max()

# Calculate RIXS at the isotropic L3 maximum. The incoming beam is linearly
# polarized; the two outgoing polarization channels are summed below.
l3_window = energy < edge_shift + 5
resonance_energy = energy[l3_window][np.argmax(xas[l3_window])]
energy_loss = np.linspace(-0.2, 10.0, 1201)
rixs = edrixs.rixs(
    eval_i, evec_i, hmat_i, hmat_n, trans_ops,
    np.array([resonance_energy]), energy_loss,
    gamma_c=gamma_c,
    gamma_f=0.05,
    thin=thin,
    thout=np.pi / 2,
    phi=phi,
    pol_type=[
        ('linear', 0, 'linear', 0),
        ('linear', 0, 'linear', np.pi / 2),
    ],
    temperature=temperature,
    backend=backend,
    backend_kws={
        'nkryl': 190,
        'linsys_tol': 1e-10,
        'linsys_maxiter': 2000,
    },
)
rixs_spectrum = rixs[0].sum(axis=-1)
rixs_spectrum /= rixs_spectrum.max()

# Plot XAS and the RIXS energy-loss spectrum side by side.
fig, (ax_xas, ax_rixs) = plt.subplots(
    1, 2, figsize=(11, 4), constrained_layout=True
)
ax_xas.plot(energy, xas)
ax_xas.axvline(resonance_energy, color='tab:red', linestyle='--', linewidth=1)
ax_xas.set_xlabel('Incident photon energy (eV)')
ax_xas.set_ylabel('Normalized intensity')
ax_xas.set_title('NiO SIAM XAS')

ax_rixs.plot(energy_loss, rixs_spectrum)
ax_rixs.set_xlabel('Energy loss (eV)')
ax_rixs.set_ylabel('Normalized intensity')
ax_rixs.set_title(r'RIXS at $E_{L_3}=%.2f$ eV' % resonance_energy)
plt.show()
