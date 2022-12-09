"""
Classes to compute arbitrary response matrices
"""
import numpy
import multiprocessing
from .elements import RMVariables, RMObservables
from .elements import sum_polab, Observable
from functools import partial
from fnmatch import fnmatch
import warnings


globring = None
globobs = None
    

class ElementResponseMatrix(object):

    def __init__(self, vkw={}, okw={}):
        self.variables = RMVariables(**vkw)
        self.observables = RMObservables(**okw)
        self.fullrm = None
        super(ElementResponseMatrix, self).__init__()
        
    def add_variables(self, ring, variables):
        self.variables.add_elements(ring, variables)
        
    def add_observables(self, ring, observables):
        self.observables.add_elements(ring, observables)
        
    def add_variables_refpts(self, ring, name, delta, refpts, attname,
                             index=None, sum_zero=False):
        self.variables.add_elements_refpts(ring, name, delta,
                                           refpts, attname, index=index)
        if sum_zero:
            assert 'Polynom' in attname,\
               'sum_zero available only for PolynomA/B attribute'
            assert index is not None,\
               'index required for sum_zero'
            sum_zero = Observable(name+'_SUM', fun=sum_polab,
                                  args=(refpts, attname, index))
            self.add_observables(ring, sum_zero)
                                                
    def add_observables_refpts(self, ring, name, refpts, index=None,
                               weight=1):
        self.observables.add_elements_refpts(ring, name, refpts,index=index,
                                             weight=weight)
                                  
    @staticmethod
    def _resp_one(ring, observables, variable):
        if ring is None:
            ring = globring
        if observables is None:
            observables = globobs
        v0 = variable.get(ring)
        variable.set(ring, v0+variable.delta)
        op = observables.values(ring)
        variable.set(ring, v0-variable.delta)  
        om = observables.values(ring)
        do = [(oom-oop)/2/variable.delta for oom, oop in zip(om, op)]
        return do         
        
    def compute_fullrm(self, ring, use_mp=False, pool_size=None, start_method=None):
        if use_mp:
            ctx = multiprocessing.get_context(start_method)
            if pool_size == None:
                pool_size = min(len(self.variables), multiprocessing.cpu_count())
            if ctx.get_start_method() == 'fork':
                global globring
                global globobs
                globring = ring
                globobs = self.observables
                args = [(None, None, var) for var in self.variables]
                with ctx.Pool(pool_size) as pool:
                    results = pool.starmap(partial(self._resp_one), args)
                globring = None
            else:
                args = [(ring, self.observables, var) for var in self.variables]
                with ctx.Pool(pool_size) as pool:
                    results = pool.starmap(partial(self._resp_one), args)
        else:
            results = [self._resp_one(ring, self.observables, var) for var in self.variables]            
        self.fullrm = numpy.array(results)         
       
    def get_mat(self):
        assert self.fullrm is not None,\
           ' Empty response matrix: please run compute_fullrm() first'                               
        return self.fullrm                          
        
    def get_vals(self, ring, masko=None, maskv=None):           
        obs = self.observables.values(ring, masko)
        var = self.variables.get(ring, maskv)
        return obs, var
        
    def _get_weights(self, mask):
        if mask is None:
            mask = numpy.ones(len(self.observables), dtype=bool)
        w = [obs.weight for obs in numpy.array(self.observables)[mask]]
        return numpy.array(w) 
        
    def svd(self, mat=None):
        if mat is None:
            mat = self.get_mat()
        return numpy.linalg.svd(mat, full_matrices=False)                                                
                
    def svd_fit(self, ring, target, svd_cut=0, apply_correction=False, masko=None, maskv=None, niter=1):
        if masko is None:
            masko = numpy.ones(len(self.observables), dtype=bool)
        if maskv is None:
            maskv = numpy.ones(len(self.variables), dtype=bool)
        mat = self.get_mat()[maskv,:][:,masko]
        weights = self._get_weights(masko)
        wm = (mat*weights).T
        dks = numpy.zeros(sum(maskv))
        for i in range(niter):
            val, var = self.get_vals(ring, masko=masko, maskv=maskv)
            err = (val-target[masko])
            U, W, V = self.svd(mat=wm)
            Winv = numpy.linalg.inv(numpy.diag(W))
            Wtmp = numpy.zeros(Winv.shape)
            Wtmp[:len(Winv)-svd_cut]=Winv[:len(Winv)-svd_cut]
            dk = V.T @ Wtmp @ U.T @ (err*weights)
            if apply_correction:
                self.variables.set(ring, var+dk, maskv)
            dks += dk
        exp = mat.T @ dks
        return dks, exp

        
class OrbitResponseMatrix(ElementResponseMatrix):

    _MAT_OBS = {'X': ['X', 'HST_SUM'],
                  'Y': ['Y', 'VST_SUM']}
    _MAT_VARS = {'X': ['HST', 'RF'], 'Y': ['VST']}         

    def __init__(self, ring, bpms, steerers, cavities=None, sum_zero=True,
                 deltax=1.0e-4, deltay=1.0e-4, **kwargs):
        super(OrbitResponseMatrix, self).__init__(**kwargs)
        self.set_bpms(ring, bpms)
        self.set_steerers(ring, steerers, deltax, deltay, sum_zero)
        if cavities is not None:
            self.set_cavities(ring, cavities)             
        
    def set_bpms(self, ring, refpts):
        self.add_observables_refpts(ring, 'closed_orbit', refpts=refpts, index=0)
        self.add_observables_refpts(ring, 'closed_orbit', refpts=refpts, index=2)
        
    def set_steerers(self, ring, refpts, deltax, deltay, sum_zero):
        self.add_variables_refpts(ring, 'HST', deltax, refpts,'PolynomB', index=0, sum_zero=sum_zero)
        self.add_variables_refpts(ring, 'VST', deltay, refpts,'PolynomA', index=0, sum_zero=sum_zero)
        
    def set_cavities(self, ring, refpts, delta=10):
        active = [e.PassMethod.endswith('CavityPass') for e in ring[refpts]]
        assert numpy.all(active), \
            'Inactive cavities used for RM, please turn on your cavities'
        self.add_variables_refpts(ring, 'RF', delta, [refpts],'Frequency')
        
    def correct(self, ring, target=0, plane=['X', 'Y'], svd_cut=0, apply_correction=False, niter=1):
        svd_cut = numpy.broadcast_to(svd_cut, len(plane))       
        target = numpy.broadcast_to(target, len(self.observables))
        dk = {}
        exp = {}
        for p, s in zip(plane, svd_cut):
            assert fnmatch(p,'[XY]'),\
               'Orbit plane has to be X or Y'
            masko = numpy.zeros(len(self.observables), dtype=bool)
            maskv = numpy.zeros(len(self.variables), dtype=bool)
            for st in self._MAT_OBS[p]:
                masko = masko | [fnmatch(obs.name, st) for obs in self.observables]
            for st in self._MAT_VARS[p]:    
                maskv = maskv | [fnmatch(var.name, st) for var in self.variables]
            dk[p], exp[p] = self.svd_fit(ring, target, svd_cut=s, apply_correction=apply_correction,
                                         masko=masko, maskv=maskv, niter=niter)
        return dk, exp 
 
        
class TrajectoryResponseMatrix(OrbitResponseMatrix): 

    def __init__(self, ring, bpms, steerers, **kwargs):
        super(TrajectoryResponseMatrix, self).__init__(ring, bpms, steerers,
                                                       sum_zero=False, okw=kwargs)    
    def set_bpms(self, ring, refpts):
        self.add_observables_refpts(ring, 'trajectory', refpts=refpts, index=0)
        self.add_observables_refpts(ring, 'trajectory', refpts=refpts, index=2)
        
    def correct(self, ring, mat=None, target=0, plane=['X', 'Y'], svd_cut=0,
                apply_correction=False, threshold=None, niter=1):
        target = numpy.broadcast_to(target, len(self.observables))
        svd_cut = numpy.broadcast_to(svd_cut, len(plane))
        threshold = numpy.broadcast_to(threshold, len(plane))
        dk = {}
        exp = {}
        for p, s, th in zip(plane, svd_cut, threshold):
            assert fnmatch(p,'[XY]'),\
               'Trajectory plane has to be X or Y'
            masko = numpy.zeros(len(self.observables), dtype=bool)
            maskv = numpy.zeros(len(self.variables), dtype=bool)
            for st in self._MAT_OBS[p]:
                masko = masko | [fnmatch(obs.name, st) for obs in self.observables]
            for st in self._MAT_VARS[p]:    
                maskv = maskv | [fnmatch(var.name, st) for var in self.variables]
            if th is not None:
                v, _ = self.get_vals(ring)
                idx_values = numpy.where(numpy.absolute(v) >= th)[0]
                if len(idx_values) > 0:
                    masko[:idx_values[0]] = False
            dk[p], exp[p] = self.svd_fit(ring, target, svd_cut=s, apply_correction=apply_correction,
                                         masko=masko, maskv=maskv, niter=niter)
        return dk, exp 
        

class OpticsResponseMatrix(ElementResponseMatrix):

    def __init__(self, ring, bpms, magnets, magnames, magattr, magidx,
                 obsnames, obsidx, obsw=1, delta=1.0e-3, fit_tune=False,
                 fit_chrom=False, **kwargs):
        super(OpticsResponseMatrix, self).__init__(okw=kwargs)
        self.set_bpms(ring, bpms, obsnames, obsidx, obsw)
        self.set_magnets(ring, magnets, magnames, magattr, magidx, delta) 
        if fit_tune:
            qx = Observable('QX', fun=self.get_tune, args=(0,), weight=len(bpms))
            qy = Observable('QY', fun=self.get_tune, args=(1,), weight=len(bpms))
            self.add_observables(ring, [qx, qy])   
        if fit_chrom:
            qpx = Observable('QPX', fun=self.get_chrom, args=(0,), weight=len(bpms))
            qpy = Observable('QPY', fun=self.get_chrom, args=(1,), weight=len(bpms))
            self.add_observables(ring, [qpx, qpy])                 
        
    @staticmethod
    def get_tune(ring, index):
        return ring.get_tune()[index]
        
    @staticmethod
    def get_chrom(ring, index):
        return ring.get_tune()[index]    
            
    def set_bpms(self, ring, refpts, obsnames, obsidx, obsw):
        assert len(obsnames) == len(obsidx), \
            'Optics RM: observables names and indexes must have the same length' 
        obsw = numpy.broadcast_to(obsw, len(obsnames))
        for o, i, w in zip(obsnames, obsidx, obsw):
            for ii in numpy.atleast_1d(i):
                self.add_observables_refpts(ring, o, refpts=refpts, index=ii, weight=w)
            
    def set_magnets(self, ring, refpts, magnames, magattr, magidx, delta):
        magnames = numpy.broadcast_to(magnames, len(refpts)) 
        magattr = numpy.broadcast_to(magattr, len(refpts))
        magidx = numpy.broadcast_to(magidx, len(refpts))
        delta = numpy.broadcast_to(delta, len(refpts))
        self.add_variables_refpts(ring, magnames, delta, refpts, magattr, index=magidx)
        
    def correct(self, ring, target, svd_cut=0, apply_correction=False, niter=1):    
        assert len(target) == len(self.observables), \
            'Linopt RM: target must have its length equal to the number of observables'
        dk, exp = self.svd_fit(ring, target, svd_cut=svd_cut, apply_correction=apply_correction, niter=1)
        return dk, exp
        
        
class LinoptResponseMatrix(OpticsResponseMatrix):

    def __init__(self, ring, bpms, quadrupoles, obsnames, obsidx, obsw=1,
                 delta=1.0e-3, fit_tune=True, **kwargs):
        super(LinoptResponseMatrix, self).__init__(ring, bpms, quadrupoles, 'QUAD',
                                                   'PolynomB', 1, obsnames, obsidx,
                                                   obsw=obsw, delta=delta,
                                                   fit_tune=fit_tune, okw=kwargs)
                                                   
                                                   
class ChromoptResponseMatrix(OpticsResponseMatrix):

    def __init__(self, ring, bpms, sextupoles, obsw=1, delta=1.0e-3,
                 fit_chrom=True, **kwargs):
        super(ChromoptResponseMatrix, self).__init__(ring, bpms, sextupoles, 'SEXT',
                                                     'PolynomB', 2, ['W'], [[0, 1]],
                                                     obsw=obsw, delta=delta,
                                                     fit_chrom=fit_chrom, okw=kwargs)                      
