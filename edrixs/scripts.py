def _run_solver(solver):
    """Call an f2py solver and detach cleanly when MPI-spawned."""
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()
    try:
        solver(fcomm, rank, size)
    finally:
        parent = MPI.Comm.Get_parent()
        if parent != MPI.COMM_NULL:
            parent.Disconnect()


def ed():
    from .fedrixs import ed_fsolver

    _run_solver(ed_fsolver)


def xas():
    from .fedrixs import xas_fsolver

    _run_solver(xas_fsolver)


def rixs():
    from .fedrixs import rixs_fsolver

    _run_solver(rixs_fsolver)


def opavg():
    from .fedrixs import opavg_fsolver

    _run_solver(opavg_fsolver)
