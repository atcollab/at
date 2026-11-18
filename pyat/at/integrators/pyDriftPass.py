"""Python implementation of DriftPass.

This is an example of a python passmethod.

The python passmethod module must be in the python path or in
*at.integrators*, and contain a function *trackFunction()*

The element creation is as follows:

drift = elements.Drift('drift', 1.0, PassMethod='pyDriftPass')

To be noted:
The particles are described by a 2-dimensional Fortran-aligned (6,num_particles) array.
In case a python passmethod and a c passmethod with the same name
are found, the python passmethod is used.
"""  # noqa: N999

import numpy as np


def drift6(r, length: float):
    p_norm = 1 / (1 + r[4])
    norm_l = length * p_norm
    r[0] += norm_l * r[1]
    r[2] += norm_l * r[3]
    r[5] += norm_l * p_norm * (r[1] * r[1] + r[3] * r[3]) / 2


# noinspection pep8-naming
def trackFunction(rin, elem=None):
    le = elem.Length
    for rtmp in np.atleast_2d(rin.T):
        if hasattr(elem, "T1"):
            rtmp[:] += elem.T1
        if hasattr(elem, "R1"):
            rtmp[:] = elem.R1 @ rtmp
        drift6(rtmp, le)
        if hasattr(elem, "R2"):
            rtmp[:] = elem.R2 @ rtmp
        if hasattr(elem, "T2"):
            rtmp[:] += elem.T2
