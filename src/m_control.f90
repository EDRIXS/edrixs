!> Runtime control parameters broadcast to all MPI ranks via config().
!!
!! All variables are declared with the save attribute so they persist across
!! subroutine calls. Default values are set in config() before reading
!! config.in.
module m_control
    use m_constants, only: dp

    implicit none

    integer, public, save :: ed_solver = 1       !< Eigensolver selection: 0=full diag, 1=Lanczos, 2=ARPACK
    integer, public, save :: num_val_orbs = 1    !< Number of valence spin-orbitals
    integer, public, save :: num_core_orbs = 1   !< Number of core spin-orbitals

    integer, public, save :: ndim_i = 1          !< Dimension of the initial Hilbert space
    integer, public, save :: ndim_n_nocore = 1   !< Dimension of intermediate space (valence sector only)
    !> Total dimension of intermediate space including core-hole degree of freedom.
    !! ndim_n = ndim_n_nocore * num_core_orbs
    integer, public, save :: ndim_n = 1
    integer, public, save :: ndim_f = 1          !< Dimension of the final Hilbert space

    integer, public, save :: nhopp_i = 1   !< Number of non-zero hopping terms for the initial Hamiltonian
    integer, public, save :: ncoul_i = 1   !< Number of non-zero Coulomb terms for the initial Hamiltonian
    integer, public, save :: nhopp_n = 1   !< Number of non-zero hopping terms for the intermediate Hamiltonian
    integer, public, save :: ncoul_n = 1   !< Number of non-zero Coulomb terms for the intermediate Hamiltonian

    integer, public, save :: ntran_xas    = 1  !< Number of non-zero XAS transition operator elements
    integer, public, save :: ntran_rixs_i = 1  !< Number of non-zero RIXS absorption transition operator elements
    integer, public, save :: ntran_rixs_f = 1  !< Number of non-zero RIXS emission transition operator elements

    integer, public, save :: neval   = 1    !< Number of lowest eigenvalues to compute
    logical, public, save :: idump   = .false.  !< Whether to write eigenvectors to disk
    integer, public, save :: nvector = 1    !< Number of eigenvectors to dump

    integer, public, save :: num_gs  = 1    !< Number of ground states used to compute XAS/RIXS

    integer, public, save :: maxiter   = 1    !< Maximum Lanczos/ARPACK iterations for the initial Hamiltonian
    integer, public, save :: linsys_max = 1  !< Maximum MINRES iterations for the linear system solver
    integer, public, save :: min_ndim  = 100 !< Minimum Hilbert-space dimension below which full diag is used
    integer, public, save :: ncv       = 1   !< Number of Arnoldi vectors for ARPACK (ncv >= neval + 2)
    integer, public, save :: nkryl     = 500 !< Dimension of the Krylov space for building XAS/RIXS spectra

    real(dp), public, save :: eigval_tol = 1E-10  !< Convergence tolerance for eigenvalues
    real(dp), public, save :: linsys_tol = 1E-12  !< Convergence tolerance for the linear system solver

    real(dp), public, save :: omega_in   !< Incident photon energy for RIXS (eV)
    real(dp), public, save :: gamma_in   !< Core-hole lifetime broadening for RIXS (eV, half-width)

    integer, public, save :: origin_nprocs = 1  !< Total number of MPI processes in the original communicator
    integer, public, save :: nprocs        = 1  !< Number of MPI processes actually used for computation
    integer, public, save :: origin_myid   = 0  !< Rank in the original communicator
    integer, public, save :: myid          = 0  !< Rank in the active communicator

    integer, public, save :: master      = 0  !< Rank of the master process
    integer, public, save :: new_comm    = 0  !< Active MPI communicator (may be a sub-communicator of origin_comm)
    integer, public, save :: origin_comm = 0  !< Original MPI communicator passed in from Python

end module m_control
