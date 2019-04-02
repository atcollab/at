"""
Load lattices from Matlab files.
"""
import scipy.io
import numpy
from . import element_from_dict


def _load_element(index, element_array, check=True, quiet=False):
    """Load what scipy produces into a pyat element object.
    """
    kwargs = {}
    for field_name in element_array.dtype.fields:
        # Remove any surplus dimensions in arrays.
        data = numpy.squeeze(element_array[field_name])
        # Convert strings in arrays back to strings.
        if data.dtype.type is numpy.unicode_:
            data = str(data)
        kwargs[field_name] = data

    return element_from_dict(kwargs, index=index, check=check, quiet=quiet)


def load_mat(filename, key=None, check=True, quiet=False):
    """Load a matlab at structure into a Python at list

    PARAMETERS
        filename        name of a '.mat' file
        key             name of the Matlab variable containing the lattice.
                        Default: Matlab variable name if there is only one,
                        otherwise 'RING'

    KEYWORDS
        check=True      if False, skip the coherence tests

    OUTPUT
        list    pyat ring
    """
    m = scipy.io.loadmat(filename)
    if key is None:
        matvars = [varname for varname in m if not varname.startswith('__')]
        key = matvars[0] if (len(matvars) == 1) else 'RING'
    element_arrays = m[key].flat
    return [_load_element(i, elem[0][0], check=check, quiet=quiet) for
            (i, elem) in enumerate(element_arrays)]
