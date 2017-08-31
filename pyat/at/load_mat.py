"""
Load lattice from Matlab file.

This is working from a specific file and may not be general.
"""

import scipy.io
import at
from at import elements
import numpy


CLASS_MAPPING = {'Quad': 'Quadrupole',
                 'Sext': 'Sextupole'}

FAMILY_MAPPING = {'AP': 'Aperture'}


def load_element(element_array):
    """
    Load what scipy produces into a pyat element object.
    """
    raw_data = element_array[0]
    kwargs = {}
    for item in element_array[0][0][0].dtype.fields:
        # Remove any surplus dimensions in arrays.
        data = numpy.squeeze(raw_data[item][0, 0])
        # Convert strings in arrays back to strings.
        if data.dtype.type is numpy.unicode_:
            data = str(data)
        kwargs[item] = data

    try:
        class_name = kwargs.pop('Class')
        class_name = CLASS_MAPPING.get(class_name, class_name)
    except KeyError:
        class_name = FAMILY_MAPPING.get(kwargs['FamName'], 'Drift')

    cl = getattr(elements, class_name)
    # Remove mandatory attributes from the keyword arguments.
    args = [kwargs.pop(attr) for attr in cl.REQUIRED_ATTRIBUTES]
    element = cl(*args, **kwargs)
    return element


def load(filename, key='RING'):
    """Load a matlab at structure into a Python at list
    """
    m = scipy.io.loadmat(filename)
    mat_ring = m[key]
    py_ring = []
    for item in mat_ring:
        py_ring.append(load_element(item))
    return py_ring


if __name__ == '__main__':
    m = load('../atmat/atmatch/ExampleATMATCH/dba.mat')
    rin = numpy.array((1e-6, 0, 0, 0, 0, 0))
    print(rin)
    at.atpass(m, rin, 1)
    print(rin)
