"""
The ctypes method.  I gave up on this because I considered it too slow.
"""
import os
from ctypes import c_double, cdll
import time

import pkg_resources
pkg_resources.require('aphla')

import numpy
import aphla as ap
ap.machines.load('SRI21')
ap.machines.use('SR')


SO_DIR = '/home/hgs15624/code/matlab/at_esrf/atintegrators'
MAX_ORDER = 3
MAX_INT_STEPS = 10


def load_functions():
    funcs = {}
    for item in os.listdir(SO_DIR):
        if item.endswith('Pass.so'):
            funcs[item[:-3]] = cdll.LoadLibrary(os.path.join(SO_DIR, item))

    for key in funcs:
        print('{0} {1}'.format(key, funcs[key].__getattr__(key)))
        f = funcs[key].__getattr__(key)
    return funcs

FUNCS = load_functions()
DRIFT_PASS_METHOD = FUNCS['DriftPass'].DriftPass
BEND_PASS_METHOD = FUNCS['BndMPoleSymplectic4E2Pass'].BndMPoleSymplectic4E2Pass
QUAD_PASS_METHOD = FUNCS['QuadLinearPass'].QuadLinearPass
SEXT_PASS_METHOD = FUNCS['StrMPoleSymplectic4Pass'].StrMPoleSymplectic4Pass


def drift_pass(rin, element):
    DRIFT_PASS_METHOD(rin.ctypes.data, c_double(element.length),
      None, None, None, None,
      None, None, 1)


def bend_pass(rin, element):
    fint1 = 0.6438
    fint2 = 0.6438
    entrance = exit = 0.0654
    bending_angle = 0.1309
    irho = bending_angle / element.length
    polyA = numpy.zeros((4,))
    polyB = numpy.zeros((4,))
    h1 = 0
    h2 = 0
    gap = 0.0466

    BEND_PASS_METHOD(rin.ctypes.data, c_double(element.length), c_double(irho),
      polyA.ctypes.data, polyB.ctypes.data,
      MAX_ORDER, MAX_INT_STEPS,
      c_double(entrance), c_double(exit),
      c_double(fint1), c_double(fint2), c_double(gap),
      c_double(h1), c_double(h2),
      None, None, None, None, 1)


def quad_pass(rin, element):
    QUAD_PASS_METHOD(rin.ctypes.data, c_double(element.length),
                     c_double(element.k1), None, None, None, None, 1)


def sext_pass(rin, element):
    polyA = numpy.zeros((4,))
    polyB = numpy.zeros((4,))
    polyB[3] = element.k2
    SEXT_PASS_METHOD(rin.ctypes.data, c_double(element.length),
                     polyA.ctypes.data, polyB.ctypes.data,
                     MAX_ORDER, MAX_INT_STEPS, None, None, None, None, 1)


function_map = {'BEND': bend_pass,
                'SEXT': sext_pass,
                'QUAD': quad_pass}


# one particle for now
def at_pass(ring, rin, num_turns):
    rout = numpy.zeros((6, num_turns))

    for i in range(num_turns):
        for element in ring:
            element.pass_method(rin, element)
            rout[:,i] = rin



if __name__ == '__main__':
    the_ring = ap.getElements('*')
    for item in the_ring:
        if item.family in function_map:
            item.pass_method = function_map[item.family]
        else:
            item.pass_method = drift_pass

    rin = numpy.zeros(6)
    rin[0] = 1e-6
    print(rin)
    t = time.time()
    at_pass(the_ring, rin, 1000)
    print('Time taken: {}'.format(time.time() - t))
    print(rin)
