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
    
_mode2range = {RMMode.PLUS : numpy.array([1]), 
               RMMode.MINUS : numpy.array([-1]),
               RMMode.PLUSMINUS : numpy.array([-1, 1])}

    
class RMVariable(Variable):
    def __init__(setfun, getfun, delta, mode=RMMode.PLUSMINUS, name=''):
        self.delta = delta*_mode2range[mode]          
        super(RMVariable, self).__init__(setfun, getfun, name=name)
        
        
class RMElementVariable(ElementVariable):
    def __init__(refpts, attrname, delta, index=None, mode=RMMode.PLUSMINUS, name=''):  
        self.delta = delta*_mode2range[mode]        
        super(RMElementVariable, self).__init__(refpts, attrname, index=index, name=name) 
        
        
class RMObservable(Constraints):
    """dummy"""


class RMElementObservable(ElementConstraints):
    """dummy"""


class RMOrbitObservable(OrbitConstraints): 
    """dummy"""


class RMLinoptObservable(LinoptConstraints):     
     """dummy"""                                        


class ResponseMatrix(object):

    def __init__(variables, observables):
        self.variables = variables
        self.observables = observables
        
    def compute(ring):
        var0 = [var.get(ring) for var in variables] 
        print(var0)
        
        
    

