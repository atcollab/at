"""
Classes to compute arbitrary response matrices
"""
from enum import Enum
import numpy
import pandas
from at.matching import Variable, ElementVariable
from at.matching import Constraints, ElementConstraints
from at.matching import OrbitConstraints, LinoptConstraints
from at.lattice import AtError, uint32_refpts


class RMMode(Enum):
    PLUS = 1
    MINUS = 2
    PLUSMINUS = 3
 
    
class ObsType(Enum)
    ORBIT = 1
    LINOPT = 2

    
class RMElementVariable(ElementVariable):
    def __init__(self, refpts, attrname, delta, index=None, name=''):  
        self.delta = delta
        self.attrname = attrname
        self.index = index
        self.refpts = refpts        
        super(RMElementVariable, self).__init__(refpts, attrname,
                                                index=index, name=name) 


class RMOrbitObservable(OrbitConstraints): 
    """dummy"""


class RMLinoptObservable(LinoptConstraints):     
     """dummy"""                                        


class ElementResponseMatrix(object):

    def __init__(self, observables, mode=RMMode.PLUSMINUS):
        self.variables = []
        self.observables = numpy.atleast_1d(observables)
        self.mode = mode
        self.fullrm = None
        self.obs_conf = None
        self.vars_conf = None
        self.excluded_refpts = None
        self.excluded_obsnames = None
        self.excluded_varnames = None
        self.excluded_attrnames = None
        self.orbitobservable = None
        self.linoptobservable = None 
        
    def add_variables(self, ring, refpts, attrname, delta, index=None):
        names = [ring[r].FamName+'_'+str(i+len(self.variables)) for i, r in enumerate(refpts)]
        self.variables += [RMElementVariable(r, attrname, delta, index=index, 
                                             name=names[i]) for i, r in enumerate(refpts)] 
                                             
    #def add_observables(self,ring, refpts, index=None, target=None, obstype=ObsType.ORBIT):
                                            
           
    
        
    def _config_vars_obs(self):
        obsnames = []
        obsref = []
        obskeys = []    
        for o in self.observables:
           if hasattr(o, 'refpts'):
               for name, refpts in zip(o.name, o.refs):               
                   for r in uint32_refpts(refpts, o.nelems):
                       obsnames.append(name)
                       obsref.append(r)
                       obskeys.append(name+'_'+str(r))
           else:
               for name in o.name:
                   obsnames.append(name)
                   obsref.append(-1)
                   obskeys.append(name)
        obsd = {}
        varsd = {} 
        dobs = {key: name for key, name in zip(obskeys, obsnames)}
        drefs = {key: ref for key, ref in zip(obskeys, obsref)}
        dattr = {v.name: v.attrname for v in self.variables}
        didx = {v.name: v.index for v in self.variables}
        obsd.update({'Obs': pandas.Series(dobs)})
        obsd.update({'Refs': pandas.Series(drefs)})
        varsd.update({'Attr': pandas.Series(dattr)})
        varsd.update({'Index': pandas.Series(didx)})
        self.obs_conf = pandas.DataFrame(obsd)
        self.vars_conf = pandas.DataFrame(varsd)                    
      
    @staticmethod             
    def flat_results(values):
        return numpy.concatenate(values, axis=1).ravel()   
        
    def set_excluded(self, refpts=None, obsnames=None, varnames=None,
                     attrnames=None):
        self.excluded_refpts = refpts
        self.excluded_obsnames = obsnames
        self.excluded_varnames = varnames
        self.excluded_attrnames = attrnames  
        
    def compute_fullrm(self,ring):
        self._config_vars_obs()
        var0 = [var.get(ring) for var in self.variables] 
        val0 = [obs.values(ring) for obs in self.observables]       
        rv = {}
        for var, v0 in zip(self.variables, var0):
            var.set(ring, v0+var.delta)  
            op = numpy.squeeze([self.flat_results(obs.values(ring))
                               for obs in self.observables])
            var.set(ring, v0-var.delta)  
            om = numpy.squeeze([self.flat_results(obs.values(ring))
                               for obs in self.observables])
            do = {key: (oom-oop)/2/var.delta for key, oom, oop 
                  in zip(self.obs_conf.index, om, op)}
            var.set(ring, v0)
            rv.update({var.name: pandas.Series(do)})
        self.fullrm = pandas.DataFrame(rv)
       
    def get_reduced_rm(self):
        rtmp = self.fullrm
        vtmp = self.vars_conf
        if self.excluded_refpts is not None:
            rtmp = rtmp.loc[~self.obs_conf.Refs.isin(self.excluded_refpts)]
        if self.excluded_obsnames is not None:
            rtmp = rtmp.loc[~self.obs_conf.Obs.isin(self.excluded_obsnames)]                                   
        if self.excluded_attrnames is not None:
            vtmp = vtmp[~vtmp.Attr.isin(self.excluded_attrnames)]              
        if self.excluded_varnames is not None:
            varnames = [v for v in vtmp.index
                        if v not in self.excluded_varnames]
        else:
            varnames = vtmp.index            
        return rtmp[varnames] 
        
    def svd(self):
        U, W, V = numpy.linalg.svd(self.get_reduced_rm().to_numpy(),
                                   full_matrices=False)
        return U, W, V
        
    def svd_fit(self, svd_cut=0):
        U, W, V = self.svd()
        Winv = numpy.linalg.inv(np.diag(self.W))
        if svd_cut>0:
            Winv[-svd_cut:] = 0
        a = numpy.dot(U.T, self.errors_vectors)
        b = numpy.dot(Winv, a)
        dk = numpy.dot(V.T, b)
        return dk       
                  
        
        
    

