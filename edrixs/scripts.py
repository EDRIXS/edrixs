def ed():
    from mpi4py import MPI
    from .fedrixs import ed_fsolver

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()

    ed_fsolver(fcomm, rank, size)


def xas():
    from mpi4py import MPI
    from .fedrixs import xas_fsolver

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()

    xas_fsolver(fcomm, rank, size)


def rixs():
    from mpi4py import MPI
    from .fedrixs import rixs_fsolver

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()

    rixs_fsolver(fcomm, rank, size)


def opavg():
    from mpi4py import MPI
    from .fedrixs import opavg_fsolver

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()

    opavg_fsolver(fcomm, rank, size)
