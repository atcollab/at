"""
Classes to compute arbitrary response matrices
"""
from enum import Enum
import numpy
import pandas
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
        self.obsnames = []
        self.buildobsnames
        
    def build_obsnames(self):
        for o in self.observables:
           if hasattr(o, 'refpts'):
               for r in o.refpts:
                   self.obsnames.append(o.name+'_'+str(r))
           else:
               self.obsnames.append(o.name)      
        
    def compute(self,ring):
        var0 = [var.get(ring) for var in self.variables] 
        val0 = [obs.values(ring) for obs in self.observables]       
             
        for var, v0 in zip(self.variables, var0):
            if self.mode == RMMode.PLUSMINUS:
                var.set(ring, v0+var.delta)  
                op = [obs.values(ring) for obs in self.observables]
                var.set(ring, v0-var.delta)  
                om = [obs.values(ring) for obs in self.observables]
                do = [(oom-oop)/2/var.delta for oom, oop in zip(om, op)]
            elif self.mode == RMMode.PLUS:
                var.set(ring, v0+var.delta)  
                op = [numpy.array(obs.values(ring)) for obs in self.observables]
                do = [(val0-oop)/var.delta for oom, oop in zip(om, op)]
            else:
                var.set(ring, v0-var.delta)  
                om = [numpy.array(obs.values(ring)) for obs in self.observables]
                do = [(oom-val0)/var.delta for oom, oop in zip(om, op)]  
            var.set(ring, v0)
            
        self.fullrm = rv
        
         
        
     

                  
        
        
    

