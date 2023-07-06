import numpy
'''
This is an example of a python passmethod, here the Drift
passmethod. This is a simplified version that does not check
aperture.

The python passmethods has to be in the python path or in 
at.integratorsand contain a function trackFunction()
'''


def trackFunction(rin, elem=None):
    scaling = elem.FieldScaling
    for i in range(0, len(rin), 6):
        rin[i+1] /= scaling
        rin[i+3] /= scaling
        rin[i+4] = (rin[i+4] + 1.0 - scaling) / scaling
    return rin
