"""
Load lattice from Matlab file.

This is working from a specific file and may not be general.
"""

import scipy.io
from .lattice import elements
import numpy


CLASS_MAPPING = {'Quad': 'Quadrupole',
                 'Sext': 'Sextupole'}

FAMILY_MAPPING = {'AP': 'Aperture'}


def load_element(element_array):
    """
    Load what scipy produces into a pyat element object.
    """
    kwargs = {}
    for field_name in element_array.dtype.fields:
        # Remove any surplus dimensions in arrays.
        data = numpy.squeeze(element_array[field_name])
        # Convert strings in arrays back to strings.
        if data.dtype.type is numpy.unicode_:
            data = str(data)
        kwargs[field_name] = data

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
    py_ring = [load_element(item[0, 0]) for item in m[key].flat]
    return py_ring
