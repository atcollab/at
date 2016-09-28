"""
The python extension method.
"""
import sys
sys.path.append('./build/lib.linux-x86_64-2.7/')
import dls_packages
import numpy
import pml
import aphla as ap
import at
import time

pml.initialise('SRI21')
the_ring = ap.getElements('*')

PASS_METHODS = {'SEXT': 'StrMPoleSymplectic4Pass',
                'QUAD': 'QuadLinearPass',
                'BEND': 'BndMPoleSymplectic4E2Pass',
                'DRIFT': 'DriftPass'}

for element in the_ring:
    family = element.family
    element.Length = element.length
    if family == 'BEND':
        element.bending_angle = 0.130899693899575
        element.entrance_angle = 0.065449846949787
        element.exit_angle = 0.065449846949787
        element.gap = 0.0466
        element.fringe_int_1 = 0.643776824034335
        element.fringe_int_2 = 0.643776824034335
        element.max_order = 3
        element.num_int_steps = 10
        element.polynom_a = numpy.array([0,0,0,0], dtype=numpy.float64)
        element.polynom_b = numpy.array([0,0,0,0], dtype=numpy.float64)
    elif family == 'SEXT':
        element.MaxOrder = 3
        element.NumIntSteps = 10
        element.PolynomA = numpy.array([0,0,0,0], dtype=numpy.float64)
        element.PolynomB = numpy.array([0,0,element.k2,0], dtype=numpy.float64)
    elif family == 'QUAD':
        element.K = element.k1
    try:
        element.PassMethod = PASS_METHODS[family]
    except KeyError:
        element.PassMethod = 'DriftPass'

rin = numpy.zeros((1, 6))
rin[0,0] = 1e-6
t = time.time()
at.atpass(the_ring, rin, 1000)
print("Time taken: {}".format(time.time() - t))
print(rin)
