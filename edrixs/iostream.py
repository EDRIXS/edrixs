__all__ = ['write_tensor', 'dump_poles', 'load_poles']

import numpy as np
import json


def write_tensor_1(tensor, fname, only_nonzeros=False, tol=1E-10, fmt_int='{:10d}',
                   fmt_float='{:.15f}'):
    (n1, ) = tensor.shape
    is_cmplx = False
    if tensor.dtype in (complex, np.complex128):
        is_cmplx = True
    space = "    "
    f = open(fname, 'w')
    for i in range(n1):
        if only_nonzeros and abs(tensor[i]) < tol:
            continue
        if is_cmplx:
            fmt_string = fmt_int + space + (fmt_float + space) * 2 + '\n'
            f.write(fmt_string.format(i + 1, tensor[i].real, tensor[i].imag))
        else:
            fmt_string = fmt_int + space + fmt_float + space + '\n'
            f.write(fmt_string.format(i + 1, tensor[i]))
    f.close()


def write_tensor_2(tensor, fname, only_nonzeros=False, tol=1E-10, fmt_int='{:10d}',
                   fmt_float='{:.15f}'):
    (n1, n2) = tensor.shape
    is_cmplx = False
    if tensor.dtype in (complex, np.complex128):
        is_cmplx = True
    space = "    "
    f = open(fname, 'w')
    for i in range(n1):
        for j in range(n2):
            if only_nonzeros and abs(tensor[i, j]) < tol:
                continue
            if is_cmplx:
                fmt_string = (fmt_int + space) * 2 + (fmt_float + space) * 2 + '\n'
                f.write(fmt_string.format(i + 1, j + 1, tensor[i, j].real, tensor[i, j].imag))
            else:
                fmt_string = (fmt_int + space) * 2 + fmt_float + space + '\n'
                f.write(fmt_string.format(i + 1, j + 1, tensor[i, j]))
    f.close()


def write_tensor_3(tensor, fname, only_nonzeros=False, tol=1E-10, fmt_int='{:10d}',
                   fmt_float='{:.15f}'):
    (n1, n2, n3) = tensor.shape
    is_cmplx = False
    if tensor.dtype in (complex, np.complex128):
        is_cmplx = True
    space = "    "
    f = open(fname, 'w')
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                if only_nonzeros and abs(tensor[i, j, k]) < tol:
                    continue
                if is_cmplx:
                    fmt_string = (fmt_int + space) * 3 + (fmt_float + space) * 2 + '\n'
                    f.write(fmt_string.format(i + 1, j + 1, k + 1, tensor[i, j, k].real,
                                              tensor[i, j, k].imag))
                else:
                    fmt_string = (fmt_int + space) * 3 + fmt_float + space + '\n'
                    f.write(fmt_string.format(i + 1, j + 1, k + 1, tensor[i, j, k]))
    f.close()


def write_tensor_4(tensor, fname, only_nonzeros=False, tol=1E-10, fmt_int='{:10d}',
                   fmt_float='{:.15f}'):
    (n1, n2, n3, n4) = tensor.shape
    is_cmplx = False
    if tensor.dtype in (complex, np.complex128):
        is_cmplx = True
    space = "    "
    f = open(fname, 'w')
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                for m in range(n4):
                    if only_nonzeros and abs(tensor[i, j, k, m]) < tol:
                        continue
                    if is_cmplx:
                        fmt_string = (fmt_int + space) * 4 + (fmt_float + space) * 2 + '\n'
                        f.write(fmt_string.format(i + 1, j + 1, k + 1, m + 1,
                                tensor[i, j, k, m].real, tensor[i, j, k, m].imag))
                    else:
                        fmt_string = (fmt_int + space) * 4 + fmt_float + space + '\n'
                        f.write(fmt_string.format(i + 1, j + 1, k + 1, m + 1, tensor[i, j, k, m]))
    f.close()


def write_tensor_5(tensor, fname, only_nonzeros=False, tol=1E-10, fmt_int='{:10d}',
                   fmt_float='{:.15f}'):
    (n1, n2, n3, n4, n5) = tensor.shape
    is_cmplx = False
    if tensor.dtype in (complex, np.complex128):
        is_cmplx = True
    space = "    "
    f = open(fname, 'w')
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                for r in range(n4):
                    for m in range(n5):
                        if only_nonzeros and abs(tensor[i, j, k, r, m]) < tol:
                            continue
                        if is_cmplx:
                            fmt_string = (fmt_int + space) * 5 + (fmt_float + space) * 2 + '\n'
                            f.write(fmt_string.format(i + 1, j + 1, k + 1, r + 1, m + 1,
                                    tensor[i, j, k, r, m].real, tensor[i, j, k, r, m].imag))
                        else:
                            fmt_string = (fmt_int + space) * 5 + fmt_float + space + '\n'
                            f.write(fmt_string.format(i + 1, j + 1, k + 1, r + 1, m + 1,
                                    tensor[i, j, k, r, m]))
    f.close()


def write_tensor(tensor, fname, only_nonzeros=False, tol=1E-10, fmt_int='{:10d}',
                 fmt_float='{:.15f}'):
    """
    Write :math:`n` -dimension numpy array to file, currently, :math:`n` can be 1, 2, 3, 4, 5.

    Parameters
    ----------
    tensor: :math:`n` d float or complex array
        The array needs to be written.
    fname: str
        File name.
    only_nonzeros: logical (default: False)
        Only write nonzero elements.
    tol: float (default: 1E-10)
        Only write the elements when their absolute value are larger than tol
        and only_nonzeros=True.
    fmt_int: str (default: '{:10d}')
        The format for printing integer numbers.
    fmt_float: str (default: '{:.15f}')
        The format for printing float numbers.

    """

    ndim = tensor.ndim
    if ndim == 1:
        write_tensor_1(tensor, fname, only_nonzeros=only_nonzeros, tol=tol,
                       fmt_int=fmt_int, fmt_float=fmt_float)
    elif ndim == 2:
        write_tensor_2(tensor, fname, only_nonzeros=only_nonzeros, tol=tol,
                       fmt_int=fmt_int, fmt_float=fmt_float)
    elif ndim == 3:
        write_tensor_3(tensor, fname, only_nonzeros=only_nonzeros, tol=tol,
                       fmt_int=fmt_int, fmt_float=fmt_float)
    elif ndim == 4:
        write_tensor_4(tensor, fname, only_nonzeros=only_nonzeros, tol=tol,
                       fmt_int=fmt_int, fmt_float=fmt_float)
    elif ndim == 5:
        write_tensor_5(tensor, fname, only_nonzeros=only_nonzeros, tol=tol,
                       fmt_int=fmt_int, fmt_float=fmt_float)
    else:
        raise Exception("error in write_tensor: ndim >5, not implemented !")


def dump_poles(obj, file_name="poles"):
    """
    Dump the objects of poles returned from XAS or RIXS calculations to file for later plotting.

    Parameters
    ----------
    obj: Python object
        Object of poles, a dict or a list of dicts.
    file_name: string
        File name.
    """
    with open(file_name+'.json', 'w') as f:
        json.dump(obj, f, indent=2)


def load_poles(file_name='poles'):
    """
    Load the objects of poles from file.

    Parameters
    ----------
    file_name: string

    Returns
    -------
    obj: Python objects
        Poles object.
    """
    with open(file_name+'.json', 'r') as f:
        obj = json.load(f)

    return obj
