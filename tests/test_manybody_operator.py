import numpy as np
import scipy.linalg
import edrixs
from edrixs.manybody_operator import (
    one_fermion_annihilation, two_fermion, four_fermion, build_opers, density_matrix
)
from edrixs.fock_basis import get_fock_bin_by_N


# --- Number operator tests ---

def test_total_number_operator_equals_N_times_identity():
    """Total number operator N̂ = Σ_i f†_i f_i equals N·I in N-electron subspace."""
    norb, noccu = 6, 2
    basis = get_fock_bin_by_N(norb, noccu)
    emat = np.eye(norb, dtype=complex)
    H = two_fermion(emat, basis)
    expected = noccu * np.eye(len(basis), dtype=complex)
    assert np.allclose(H, expected)


def test_total_number_operator_for_three_electrons():
    """Number operator correctly gives N=3 for 3-electron basis."""
    norb, noccu = 6, 3
    basis = get_fock_bin_by_N(norb, noccu)
    emat = np.eye(norb, dtype=complex)
    H = two_fermion(emat, basis)
    assert np.allclose(H, noccu * np.eye(len(basis), dtype=complex))


# --- Hermiticity tests ---

def test_two_fermion_hermitian_when_emat_hermitian():
    """two_fermion produces Hermitian H when emat is Hermitian."""
    norb, noccu = 4, 2
    basis = get_fock_bin_by_N(norb, noccu)
    np.random.seed(0)
    A = np.random.random((norb, norb)) + 1j * np.random.random((norb, norb))
    emat = A + np.conj(A.T)
    H = two_fermion(emat, basis)
    assert np.allclose(H, np.conj(H.T))


def test_four_fermion_hermitian_from_slater_tensor():
    """four_fermion produces Hermitian H when umat comes from Slater integrals."""
    norb, noccu = 6, 2
    basis = get_fock_bin_by_N(norb, noccu)
    umat = edrixs.get_umat_slater('t2g', 3.0, 1.5, 0.5)
    H = four_fermion(umat, basis)
    assert np.allclose(H, np.conj(H.T))


def test_build_opers_2fermion_hermitian():
    """build_opers with nfermion=2 and Hermitian emat gives Hermitian matrix."""
    norb, noccu = 4, 2
    basis = get_fock_bin_by_N(norb, noccu)
    np.random.seed(1)
    A = np.random.random((norb, norb)) + 1j * np.random.random((norb, norb))
    emat = A + np.conj(A.T)
    H = build_opers(2, emat, basis)
    assert np.allclose(H, np.conj(H.T))


def test_build_opers_4fermion_hermitian():
    """build_opers with nfermion=4 and Slater tensor gives Hermitian matrix."""
    norb, noccu = 6, 2
    basis = get_fock_bin_by_N(norb, noccu)
    umat = edrixs.get_umat_slater('t2g', 3.0, 1.5, 0.5)
    H = build_opers(4, umat, basis)
    assert np.allclose(H, np.conj(H.T))


# --- Zero coefficient tests ---

def test_two_fermion_zero_emat_gives_zero_hamiltonian():
    """two_fermion with zero emat gives zero Hamiltonian."""
    norb, noccu = 4, 2
    basis = get_fock_bin_by_N(norb, noccu)
    emat = np.zeros((norb, norb), dtype=complex)
    H = two_fermion(emat, basis)
    assert np.allclose(H, 0)


def test_four_fermion_zero_umat_gives_zero_hamiltonian():
    """four_fermion with zero umat gives zero Hamiltonian."""
    norb, noccu = 4, 2
    basis = get_fock_bin_by_N(norb, noccu)
    umat = np.zeros((norb, norb, norb, norb), dtype=complex)
    H = four_fermion(umat, basis)
    assert np.allclose(H, 0)


# --- Shape tests ---

def test_two_fermion_output_shape():
    """two_fermion output has shape (len_basis, len_basis)."""
    norb, noccu = 6, 2
    basis = get_fock_bin_by_N(norb, noccu)
    emat = np.eye(norb, dtype=complex)
    H = two_fermion(emat, basis)
    n = len(basis)
    assert H.shape == (n, n)


def test_four_fermion_output_shape():
    """four_fermion output has shape (len_basis, len_basis)."""
    norb, noccu = 6, 2
    basis = get_fock_bin_by_N(norb, noccu)
    umat = edrixs.get_umat_slater('t2g', 3.0, 1.5, 0.5)
    H = four_fermion(umat, basis)
    n = len(basis)
    assert H.shape == (n, n)


def test_build_opers_shape_matches_basis():
    """build_opers output shape matches Fock basis size."""
    norb, noccu = 4, 2
    basis = get_fock_bin_by_N(norb, noccu)
    emat = np.eye(norb, dtype=complex)
    H = build_opers(2, emat, basis)
    n = len(basis)
    assert H.shape == (n, n)


def test_build_opers_batch_2fermion_shape():
    """build_opers with batch of emat returns correct shape."""
    norb, noccu = 4, 2
    basis = get_fock_bin_by_N(norb, noccu)
    emat_batch = np.eye(norb, dtype=complex)[np.newaxis, :, :].repeat(3, axis=0)
    H = build_opers(2, emat_batch, basis)
    n = len(basis)
    assert H.shape == (3, n, n)


# --- Eigenvalue tests ---

def test_full_hamiltonian_eigenvalues_real():
    """Eigenvalues of the full many-body Hamiltonian are real (Hermitian)."""
    norb, noccu = 6, 2
    basis = get_fock_bin_by_N(norb, noccu)
    cf = edrixs.cf_cubic_d(1.0)
    umat = edrixs.get_umat_slater('t2g', 3.0, 1.5, 0.5)
    H = build_opers(2, cf[0:6, 0:6], basis) + build_opers(4, umat, basis)
    evals = np.linalg.eigvalsh(H)
    assert np.allclose(evals.imag, 0)


def test_ground_state_energy_is_lowest():
    """Ground state energy is less than or equal to all other eigenvalues."""
    norb, noccu = 6, 2
    basis = get_fock_bin_by_N(norb, noccu)
    cf = edrixs.cf_cubic_d(2.0)
    H = build_opers(2, cf[0:6, 0:6], basis)
    evals = scipy.linalg.eigh(H, eigvals_only=True)
    assert evals[0] <= evals[-1]


# --- density_matrix tests ---

def test_density_matrix_shape():
    """density_matrix output has shape (len_basis, len_basis)."""
    norb, noccu = 4, 2
    basis = get_fock_bin_by_N(norb, noccu)
    dm = density_matrix(0, 0, basis, basis)
    n = len(basis)
    assert dm.shape == (n, n)


def test_density_matrix_diagonal_is_occupation():
    """Diagonal of density_matrix gives orbital occupation numbers."""
    norb, noccu = 4, 2
    basis = get_fock_bin_by_N(norb, noccu)
    # Orbital 0: count how many basis states have orbital 0 occupied
    dm = density_matrix(0, 0, basis, basis)
    expected_trace = sum(state[0] for state in basis)
    assert np.isclose(np.trace(dm).real, expected_trace)


# --- one_fermion_annihilation tests ---

def test_one_fermion_annihilation_shape():
    """one_fermion_annihilation output has correct shape."""
    lb = get_fock_bin_by_N(4, 1)  # N-1 particle basis
    rb = get_fock_bin_by_N(4, 2)  # N particle basis
    ann = one_fermion_annihilation(0, lb, rb)
    assert ann.shape == (len(lb), len(rb))


def test_one_fermion_annihilation_removes_one_particle():
    """Annihilation operator reduces particle number by 1."""
    norb, noccu = 4, 2
    basis_N = get_fock_bin_by_N(norb, noccu)
    basis_Nm1 = get_fock_bin_by_N(norb, noccu - 1)
    # Each column should connect a 2-particle state to a 1-particle state
    ann = one_fermion_annihilation(0, basis_Nm1, basis_N)
    # Matrix should be non-zero only in columns where orbital 0 is occupied
    for j, state in enumerate(basis_N):
        if state[0] == 0:
            assert np.allclose(ann[:, j], 0)
        else:
            assert not np.allclose(ann[:, j], 0)


# --- Independent regression oracles ---

def test_one_fermion_annihilation_exact_matrix_and_fermionic_sign():
    """Pin both destination rows and signs for an interior orbital."""
    lb = get_fock_bin_by_N(3, 1)
    rb = get_fock_bin_by_N(3, 2)
    expected = np.array([
        [-1, 0, 0],
        [0, 0, 0],
        [0, 0, 1],
    ], dtype=complex)
    assert np.allclose(one_fermion_annihilation(1, lb, rb), expected)


def test_two_fermion_single_hop_exact_matrix_element():
    """A single coefficient must land in the correct many-body matrix entry."""
    basis = get_fock_bin_by_N(3, 1)
    emat = np.zeros((3, 3), dtype=complex)
    emat[0, 1] = 2.0 + 0.5j
    expected = np.zeros((3, 3), dtype=complex)
    expected[0, 1] = emat[0, 1]
    assert np.allclose(two_fermion(emat, basis), expected)


def test_density_matrix_offdiagonal_exact_sign_and_destination():
    """Pin an off-diagonal density operator where the fermionic sign is negative."""
    basis = get_fock_bin_by_N(3, 2)
    expected = np.zeros((3, 3), dtype=complex)
    expected[0, 2] = -1.0
    assert np.allclose(density_matrix(0, 2, basis, basis), expected)


def test_four_fermion_exact_single_matrix_element():
    """A nonzero four-fermion tensor term must not collapse to a zero matrix."""
    basis = get_fock_bin_by_N(4, 2)
    umat = np.zeros((4, 4, 4, 4), dtype=complex)
    umat[2, 3, 1, 0] = 1.75
    expected = np.zeros((len(basis), len(basis)), dtype=complex)
    expected[-1, 0] = 1.75
    assert np.allclose(four_fermion(umat, basis), expected)


def test_known_diagonal_one_body_spectrum():
    """Use a model with an analytically known two-particle spectrum."""
    basis = get_fock_bin_by_N(4, 2)
    emat = np.diag([-3.0, -1.0, 2.0, 7.0])
    evals = np.sort(np.linalg.eigvalsh(two_fermion(emat, basis)))
    expected = np.sort([-4.0, -1.0, 4.0, 1.0, 6.0, 9.0])
    assert np.allclose(evals, expected)
