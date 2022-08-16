"""
Classes to compute arbitrary response matrices
"""
from enum import Enum
import numpy
import pandas
#  needed -> ElementVariable.setv is not pickable
import multiprocess as multiprocessing
from functools import partial
from at.matching import Variable, ElementVariable
from at.matching import Constraints, ElementConstraints
from at.matching import OrbitConstraints, LinoptConstraints
from at.lattice import AtError, uint32_refpts


_ORBIT_NAMES = ['X', 'XP', 'Y', 'YP', 'DP', 'CT']

globring = None

            
def _flat_results(values):
    return numpy.concatenate(values, axis=1).ravel() 

   
class ObsType(Enum):
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

class RMObservable(Constraints): 
    """dummy""" 


class RMOrbitObservable(OrbitConstraints): 
    """dummy"""


class RMLinoptObservable(LinoptConstraints):     
     """dummy"""                                        


class ElementResponseMatrix(object):

    def __init__(self, ring, dp=0.0):
        self.ring = ring.copy()
        self.dp =dp
        #  add attribute to ring such that they are
        #  accessible through the global variable
        self.ring.variables = []
        self.ring.observables = []
        self.ring.obs_conf = None
        self.ring.vars_conf = None
        self.ring.arbobservable = RMObservable()
        self.ring.orbitobservable = RMOrbitObservable(ring, dp) 
        self.ring.linoptobservable = RMLinoptObservable(ring, dp=0.0)
        self.fullrm = None
        self.excluded_refpts = None
        self.excluded_obsnames = None
        self.excluded_varnames = None
        self.excluded_attrnames = None
        self.excluded_obskeys = None  
        
    def add_variables(variables):
        variables = numpy.atleast_1d(variables)
        for v in variables:
            v.name = v.name+'.'+str(len(self.ring.variables)+1)
            self.ring.variables += [v]
            
    def add_observable(fun, target, name, weight=1):
        self.ring.arbobservables.add(fun, target, name=name, weight=1)
        
    def add_element_variables(self, refpts, attrname, delta, index=None, name=None):
        refpts = numpy.atleast_1d(refpts)
        if name is None:
            names = [self.ring[r].FamName+'.'+str(i+len(self.ring.variables))
                     for i, r in enumerate(refpts)]
        else:
            names = numpy.broadcast_to(name, (len(refpts), ))
            names = [names[i]+'.'+str(i+len(self.ring.variables))
                     for i in range(len(refpts))]
        self.ring.variables += [RMElementVariable(r, attrname, delta, index=index, 
                                                  name=names[i])
                                for i, r in enumerate(refpts)] 
                                             
    def add_element_observables(self, refpts, index, target=None, obstype=ObsType.ORBIT):
        if obstype is ObsType.ORBIT:
            self._add_orbit_observables(refpts, index, target=target)
        elif obstype is ObsType.LINOPT:
            self._add_linopt_observables(refpts, index, target=target)
        else:
            raise AtError('ObsType {} not defined.'.format(obstype))
            
    def _add_orbit_observables(self, refpts, index, target=None): 
        for r in refpts:         
            self.ring.orbitobservable.add(target, refpts=r, index=[index],
                                          name=_ORBIT_NAMES[index])           
    
    def _add_linopt_observables(self, refpts, index, target=None):
        """Empty"""
        
    def _fillobservables(self):
        if len(self.ring.linoptobservable.name)>0:
            self.ring.observables.append(self.ring.linoptobservable)
        if len(self.ring.orbitobservable.name)>0:
            self.ring.observables.append(self.ring.orbitobservable)  
        if len(self.ring.arbobservable.name)>0:
            self.ring.observables.append(self.ring.arbobservable) 
        
    def _config_vars_obs(self):
        self._fillobservables()
        obsnames = []
        obsref = []
        obskeys = []  
        obstarget = [] 
        obsvals = []   
        for o in self.ring.observables:
           vals = o.values(self.ring)              
           for name, refpts, target, v in zip(o.name, o.refs, o.target, vals):  
               if numpy.any(refpts):                
                   for r, t in zip(uint32_refpts(refpts, o.nelems), target):
                       obsnames.append(name)
                       obsref.append(r)
                       obstarget.append(t)
                       kname = self.ring[r].FamName+'.'+name+'.'+str(r)
                       obskeys.append(kname)
                       obsvals.append(numpy.squeeze(v))
               else:
                   obsnames.append(name)
                   obstarget.append(target)
                   obskeys.append(name)
                   obsvals.append(numpy.squeeze(v))
               
        obsd = {}
        varsd = {} 
        dobs = {key: name for key, name in zip(obskeys, obsnames)}
        drefs = {key: ref for key, ref in zip(obskeys, obsref)}
        dtarget = {key: target for key, target in zip(obskeys, obstarget)}
        dvals = {key: v for key, v in zip(obskeys, obsvals)}
        dattr = {v.name: v.attrname for v in self.ring.variables}
        didx = {v.name: v.index for v in self.ring.variables}
        obsd.update({'Obs': pandas.Series(dobs)})
        obsd.update({'Refs': pandas.Series(drefs)})
        obsd.update({'Target': pandas.Series(dtarget)})
        obsd.update({'Value': pandas.Series(dvals)})
        varsd.update({'Attr': pandas.Series(dattr)})
        varsd.update({'Index': pandas.Series(didx)})
        self.ring.obs_conf = pandas.DataFrame(obsd)
        self.ring.vars_conf = pandas.DataFrame(varsd)                      
        
    def set_excluded(self, refpts=None, obskeys=None, obsnames=None,
                     varnames=None, attrnames=None):
        if refpts is not None:
            self.excluded_refpts = refpts
        if obskeys is not None:
            self.excluded_obskeys = obskeys
        if obsnames is not None:
            self.excluded_obsnames = obsnames
        if varnames is not None:
            self.excluded_varnames = varnames
        if attrnames is not None:
            self.excluded_attrnames = attrnames        
        
    def set_target(obslist, targetlist):
        for o, t in zip(obslist, targetlist):
            if isinstance(o, int):
                self.ring.obs_conf.iloc[o].Target = t
            else:
                self.ring.obs_conf.loc[o].Target = t  
                
    @staticmethod
    def _resp_one(ring, variable):
        if ring is None:
            ring = globring
        v0 = variable.get(ring)
        variable.set(ring, v0+variable.delta)  
        op = numpy.squeeze([_flat_results(obs.values(ring))
                            for obs in ring.observables])
        variable.set(ring, v0-variable.delta)  
        om = numpy.squeeze([_flat_results(obs.values(ring))
                            for obs in ring.observables])
        do = {key: (oom-oop)/2/variable.delta for key, oom, oop 
              in zip(ring.obs_conf.index, om, op)}
        variable.set(ring, v0)
        return {variable.name: pandas.Series(do)}         
        
    def compute_fullrm(self, use_mp=False, pool_size=None, start_method=None):
        self._config_vars_obs()
        rv = {}
        if use_mp:
            ctx = multiprocessing.get_context(start_method)
            if pool_size == None:
                pool_size = min(len(self.ring.variables),
                                multiprocessing.cpu_count())
            if ctx.get_start_method() == 'fork':
                global globring
                globring = self.ring
                args = [(None, var) for var in self.ring.variables]
                with ctx.Pool(pool_size) as pool:
                    results = pool.starmap(partial(self._resp_one), args)
                globring = None
            else:
                args = [(self.ring, var) for var in self.ring.variables]
                with ctx.Pool(pool_size) as pool:
                    results = pool.starmap(partial(self._resp_one), args)
            for r in results:
                rv.update(r)             
        else: 
            for var in self.ring.variables:  
                rv.update(self.resp_one(self.ring, var,
                                        self.ring.observables,
                                        self.ring.obs_conf))
        self.fullrm = pandas.DataFrame(rv)
       
    def get_reduced_rm(self):
        rtmp = self.fullrm
        vtmp = self.ring.vars_conf
        if self.excluded_obskeys is not None:
            rtmp = rtmp.loc[[r for r in rtmp.index if r not in self.excluded_obskeys]]
        if self.excluded_refpts is not None:
            rtmp = rtmp.loc[~self.ring.obs_conf.Refs.isin(self.excluded_refpts)]
        if self.excluded_obsnames is not None:
            rtmp = rtmp.loc[~self.ring.obs_conf.Obs.isin(self.excluded_obsnames)]                                   
        if self.excluded_attrnames is not None:
            vtmp = vtmp[~vtmp.Attr.isin(self.excluded_attrnames)]
        varnames = vtmp.index              
        if self.excluded_varnames is not None:
            varnames = [v for v in vtmp.index if v not in self.excluded_varnames]                        
        return rtmp[varnames]
        
    def get_reduced_vals(self):
        obstmp = self.ring.obs_conf
        if self.excluded_obskeys is not None:
            obstmp = obstmp.loc[[o for o in obstmp.index if o not in self.excluded_obskeys]]
        if self.excluded_refpts is not None:
            obstmp = obstmp.loc[~self.ring.obs_conf.Refs.isin(self.excluded_refpts)]
        if self.excluded_obsnames is not None:
            obstmp = obstmp.loc[~self.ring.obs_conf.Obs.isin(self.excluded_obsnames)]                                               
        return obstmp.Value, obstmp.Target
        
    def svd_inv(self):
        mat = self.get_reduced_rm()
        self.U, self.W, self.V = numpy.linalg.svd(mat, full_matrices=False)
        self.Winv = numpy.linalg.inv(numpy.diag(self.W))      
                
    def svd_fit(self, svd_cut=0):
        Wtmp = numpy.zeros(self.Winv.shape)
        Wtmp[:len(self.Winv)-svd_cut]=self.Winv[:len(self.Winv)-svd_cut]
        val, target = self.get_reduced_vals()
        err = (val-target)
        return self.V.T @ Wtmp @ self.U.T @ err
        
    def apply_corrections(self, dk):
        var0 = [var.get(self.ring) for var in self.variables]
        for var, v0 in zip(self.variables, var0, dk):
            var.set(self.ring, v0-dk)
                    
    def get_fitted_obs(self, dk):
        mat = self.get_reduced_rm()
        return mat @ dk
