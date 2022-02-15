import numpy

'''
This is an example of a python passmethod, here the Identity
passmethod.

The python passmethods has to be in the python path or in
at.integratorsand contain a function trackFunction()

The element creation is as follows:

params = {'Length':0,
          'PassMethod':'pyIdentityPass',
          }

ib_elem = Element('py_id',**params)

To be noted:
The particles are described by a 2 dimensional vector
(6,num_particles)
In case a pyhton passmethod and a c passmethod with the same name
are found, the c passmethod is used
'''


_id = numpy.identity(6)


def trackFunction(rin, elem=None):
    for rtmp in numpy.atleast_2d(rin.T):
        rtmp[:] = numpy.dot(_id, rtmp)
