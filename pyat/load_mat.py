"""
Load lattice from Matlab file.

This is working from a specific file and may not be general.
"""

import scipy.io
import pyat
from pyat import elements
import numpy


SCALAR_FIELDS = ['Length', 'K', 'BendingAngle', 'EntranceAngle', 'ExitAngle',
                 'MaxOrder', 'NumIntSteps', 'Voltage', 'Frequency', 'HarmNumber', 'TimeLag', 'Energy']


def extract_scalars(kwargs):
    for item in kwargs:
        if item in SCALAR_FIELDS:
            kwargs[item] = kwargs[item][0]
    return kwargs


def load_element(element_array):
    """
    Load what scipy produces into a pyat element object.
    """
    data = element_array[0]
    kwargs = {}
    for item in element_array[0][0][0].dtype.fields:
        kwargs[item] = data[item][0, 0][0]

    kwargs = extract_scalars(kwargs)
    class_name = kwargs.pop('Class')
    cl = getattr(elements, class_name)
    # Remove mandatory attributes from the keyword arguments.
    args = [kwargs.pop(attr) for attr in cl.REQUIRED_ATTRIBUTES]
    element = cl(*args, **kwargs)
    return element


def load(filename):
    m = scipy.io.loadmat(filename)
    mat_ring = m['RING']
    py_ring = []
    for item in mat_ring:
        py_ring.append(load_element(item))
    return py_ring


if __name__ == '__main__':
    m = load('../atmat/atmatch/ExampleATMATCH/dba.mat')
    rin = numpy.array((1e-6, 0, 0, 0, 0, 0))
    print(rin)
    pyat.atpass(m, rin, 1)
    print(rin)
