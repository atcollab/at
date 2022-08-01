"""
Classes to compute arbitrary response matrices
"""
from enum import Enum
import numpy
from at.matching import Variable, ElementVariable
from at.matching import Constraints, ElementConstraints
from at.matching import OrbitConstraints, LinoptConstraints
from at.lattice import AtError


class RMMode(Enum):
    PLUS = 1
    MINUS = 2
    PLUSMINUS = 3

    
class RMVariable(Variable):
    def __init__(self, setfun, getfun, delta, name=''):
        self.delta = delta          
        super(RMVariable, self).__init__(setfun, getfun,
                                         name=name)
        
        
class RMElementVariable(ElementVariable):
    def __init__(self, refpts, attrname, delta, index=None, name=''):  
        self.delta = delta        
        super(RMElementVariable, self).__init__(refpts, attrname,
                                                index=index, name=name) 
        
        
class RMObservable(Constraints):
    """dummy"""


class RMElementObservable(ElementConstraints):
    """dummy"""


class RMOrbitObservable(OrbitConstraints): 
    """dummy"""


class RMLinoptObservable(LinoptConstraints):     
     """dummy"""                                        


class ResponseMatrix(object):

    def __init__(self, variables, observables, mode=RMMode.PLUSMINUS):
        self.variables = numpy.atleast_1d(variables)
        self.observables = numpy.atleast_1d(observables)
        self.mode = mode
        self.fullrm = None
        
    def compute(self,ring):
        var0 = [var.get(ring) for var in self.variables] 
        val0 = [obs.values(ring) for obs in self.observables]       
        rv = []     
        for var, v0 in zip(self.variables, var0):
            if self.mode == RMMode.PLUSMINUS:
                var.set(ring, v0 + var.delta)  
                op = numpy.squeeze([obs.values(ring) for obs in self.observables])
                var.set(ring, v0 - var.delta)  
                om = numpy.squeeze([obs.values(ring) for obs in self.observables])
                rv.append((om-op)/2/var.delta)
            elif self.mode == RMMode.PLUS:
                var.set(ring, v0 + var.delta)  
                op = numpy.squeeze([obs.values(ring) for obs in self.observables])
                rv.append((val0-op)/var.delta)
            else:
                var.set(ring, v0 - var.delta)  
                om = numpy.squeeze([obs.values(ring) for obs in self.observables])
                rv.append((om-val0)/var.delta)    
            var.set(ring, v0)
        self.fullrm = numpy.array(rv)
        
     def get_reduced_rm(self, mask_var, mask_obs, mask_refpts):
        hasrefpts = [hasattr(obs, 'refpts'
         
        
     

                  
        
        
    

