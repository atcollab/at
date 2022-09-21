"""
Classes to compute arbitrary response matrices
"""
import numpy
import pandas
#  needed -> ElementVariable.setv is not pickable
import multiprocess as multiprocessing
from .elements import RMVariables, RMObservables
from .elements import sum_polab, Observable, load_confs
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
        self.excl_obs = []
        self.excl_var = []
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
    def _resp_one(ring, observables, variable, key):
        if ring is None:
            ring = globring
        if observables is None:
            observables = globobs
        v0 = variable.get(ring)
        variable.set(ring, v0+variable.delta)
        op = observables.values(ring)
        variable.set(ring, v0-variable.delta)  
        om = observables.values(ring)
        do = {key: (oom-oop)/2/variable.delta for key, oom, oop 
              in zip(observables.conf.index, om, op)}
        variable.set(ring, v0)
        return {key: pandas.Series(do)}         
        
    def compute_fullrm(self, ring, use_mp=False, pool_size=None, start_method=None):
        rv = {}
        if use_mp:
            ctx = multiprocessing.get_context(start_method)
            if pool_size == None:
                pool_size = min(len(self.variables), multiprocessing.cpu_count())
            if ctx.get_start_method() == 'fork':
                global globring
                global globobs
                globring = ring
                globobs = self.observables
                args = [(None, None, var, key) 
                        for var, key in zip(self.variables, self.variables.conf.index)]
                with ctx.Pool(pool_size) as pool:
                    results = pool.starmap(partial(self._resp_one), args)
                globring = None
            else:
                args = [(ring, self.observables, var, key)
                        for var, key in zip(self.variables, self.variables.conf.index)]
                with ctx.Pool(pool_size) as pool:
                    results = pool.starmap(partial(self._resp_one), args)
            for r in results:
                rv.update(r)             
        else: 
            for var, key in zip(self.variables, self.variables.conf.index):  
                rv.update(self._resp_one(ring, self.observables, var, key))
        self.fullrm = pandas.DataFrame(rv)         
        
    def save_fullrm(self, filename):
        assert self.fullrm is not None,\
           ' Empty response matrix: please run compute_fullrm() first' 
        if not filename.endswith('.h5'):
            filename += '.h5'
        store = pandas.HDFStore(filename)
        store['fullrm'] = self.fullrm
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', pandas.errors.PerformanceWarning)
            store['vconf'] = self.variables.conf
            store['oconf'] = self.observables.conf
        store.close()
        
    def load_fullrm(self, ring, filename):
        store = pandas.HDFStore(filename)
        self.fullrm = store['fullrm']
        vconf = store['vconf']
        oconf = store['oconf']
        store.close()
        load_confs(self, ring, oconf, vconf)
        
    def exclude(self, obsdictlist=None, vardictlist=None):
    
        def remove_element(conf, df):
            cond = False
            for c in numpy.atleast_1d(conf):
                cond2 =True
                for k, v in c.items():
                    cond2 = cond2 & numpy.array([d[0] in v for d in df[k]], dtype=bool)
            cond = cond | cond2
            return df.loc[cond].index   
        
        if vardictlist is not None:           
            self.excl_var = remove_element(vardictlist, self.variables.conf)
        if obsdictlist is not None:
            self.excl_obs = remove_element(obsdictlist, self.observables.conf)
       
    def get_mat(self):
        assert self.fullrm is not None,\
           ' Empty response matrix: please run compute_fullrm() first'                               
        return self.fullrm.loc[~self.fullrm.index.isin(self.excl_obs),
                               ~self.fullrm.columns.isin(self.excl_var)]                           
        
    def get_vals(self, ring, mat=None):
        if mat is None:
            mat = self.get_mat()
        mask = self.observables.conf.index.isin(mat.index)  
        obs = numpy.array(self.observables.values(ring, mask))
        mask = self.variables.conf.index.isin(mat.columns)
        var = numpy.array(self.variables.get(ring, mask))
        return obs, var
        
    def svd(self, mat=None):
        if mat is None:
            mat = self.get_mat()
        return numpy.linalg.svd(mat, full_matrices=False)                                                
                
    def svd_fit(self, ring, target, mat=None, svd_cut=0, apply_correction=False):
        if mat is None:
            mat = self.get_mat()
        val, var = self.get_vals(ring, mat=mat)
        err = (val-target)
        U, W, V = self.svd(mat=mat)
        Winv = numpy.linalg.inv(numpy.diag(W))
        Wtmp = numpy.zeros(Winv.shape)
        Wtmp[:len(Winv)-svd_cut]=Winv[:len(Winv)-svd_cut]
        dk = V.T @ Wtmp @ U.T @ err
        exp = mat @ dk
        corr = {'Vals': pandas.Series(var, mat.columns),
                'Corr': pandas.Series(dk, mat.columns),
                'Exp': pandas.Series(var+dk, mat.columns)}
        obs = {'Vals': pandas.Series(val, mat.index),
               'Target': pandas.Series(target, mat.index),
               'Err': pandas.Series(err, mat.index),
               'Exp': pandas.Series(exp, mat.index)}
        if apply_correction:
            mask = self.variables.conf.index.isin(mat.columns) 
            self.variables.set(ring, var+dk, mask)
        valn, varn = self.get_vals(ring, mat=mat) 
        corr.update({'New': pandas.Series(varn, mat.columns)})
        obs.update({'New': pandas.Series(valn, mat.index)})
        return pandas.DataFrame(corr), pandas.DataFrame(obs)

        
class OrbitResponseMatrix(ElementResponseMatrix):

    _MAT_INDEX = {'X': ['X', 'HST_SUM'],
                  'Y': ['Y', 'VST_SUM']}
    _MAT_COLS = {'X': ['HST', 'RF'], 'Y': ['VST']}         

    def __init__(self, ring, bpms, steerers, cavities=None, sum_zero=True,
                 deltax=1.0e-4, deltay=1.0e-4, **kwargs):
        super(OrbitResponseMatrix, self).__init__(**kwargs)
        self.set_bpms(ring, bpms)
        self.set_steerers(ring, steerers, deltax=deltax, deltay=deltay, sum_zero=sum_zero)
        if cavities is not None:
            self.set_cavities(ring, cavities)             
        
    def set_bpms(self, ring, refpts):
        self.add_observables_refpts(ring, 'closed_orbit', refpts=refpts, index=0)
        self.add_observables_refpts(ring, 'closed_orbit', refpts=refpts, index=2)
        
    def set_steerers(self, ring, refpts, deltax=1.0e-4, deltay=1.0e-4, sum_zero=True):
        self.add_variables_refpts(ring, 'HST', deltax, refpts,'PolynomB', index=0, sum_zero=sum_zero)
        self.add_variables_refpts(ring, 'VST', deltay, refpts,'PolynomA', index=0, sum_zero=sum_zero)
        
    def set_cavities(self, ring, refpts, delta=10):
        active = [e.PassMethod.endswith('CavityPass') for e in ring[refpts]]
        assert numpy.all(active), \
            'Inactive cavities used for RM, please turn on your cavities'
        self.add_variables_refpts(ring, 'RF', delta, [refpts],'Frequency')
        
    def set_excluded_bpms(self, bpmdictlist):
        self.set_exluded(obsdistlist=bpmdictlist)
        
    def set_excluded_steerer(self, steerersdictlist):
        self.set_exluded(vardistlist=steerersdictlist)
        
    def correct(self, ring, mat=None, target=0, plane=['X', 'Y'], svd_cut=0, apply_correction=False):
        if mat is None:
            mat = self.get_mat()
        svd_cut = numpy.broadcast_to(svd_cut, len(plane))       
        target = numpy.broadcast_to(target, mat.shape[0])
        corr = pandas.DataFrame()
        obs = pandas.DataFrame()
        for p, s in zip(plane, svd_cut):
            assert fnmatch(p,'[XY]'),\
               'Orbit plane has to be X or Y'
            mati = mat.loc[mat.index.isin(self._MAT_INDEX[p], level=0),
                           mat.columns.isin(self._MAT_COLS[p], level=0)]
            ti = target[mat.index.isin(self._MAT_INDEX[p], level=0)]
            c, o = self.svd_fit(ring, ti, mat=mati, svd_cut=s, apply_correction=apply_correction)
            corr = pandas.concat((corr, c))
            obs = pandas.concat((obs, o))
        return corr, obs
 
        
class TrajectoryResponseMatrix(OrbitResponseMatrix): 

    def __init__(self, ring, bpms, steerers, **kwargs):
        super(TrajectoryResponseMatrix, self).__init__(ring, bpms, steerers,
                                                       sum_zero=False, okw=kwargs)    
    def set_bpms(self, ring, refpts):
        self.add_observables_refpts(ring, 'trajectory', refpts=refpts, index=0)
        self.add_observables_refpts(ring, 'trajectory', refpts=refpts, index=2)
        
    def correct(self, ring, mat=None, target=0, plane=['X', 'Y'], svd_cut=0,
                           apply_correction=False, threshold=None):
        if mat is None:
            mat = self.get_mat()
        target = numpy.broadcast_to(target, mat.shape[0])
        svd_cut = numpy.broadcast_to(svd_cut, len(plane))
        threshold = numpy.broadcast_to(threshold, len(plane))
        corr = pandas.DataFrame()
        obs = pandas.DataFrame()
        for p, s, th in zip(plane, svd_cut, threshold):
            assert fnmatch(p,'[XY]'),\
               'Trajectory plane has to be X or Y'
            mask = [fnmatch(idx, p+'*') for idx in mat.index.get_level_values(0)]
            mati = mat.loc[mask, mat.columns.isin(self._MAT_COLS[p], level=0)]
            ti = target[mask]
            if th is not None:
                v, _ = self.get_vals(ring, mat=mati)
                idx_values = numpy.where(numpy.absolute(v) >= th)[0]
                if len(idx_values) > 0:
                    mati = mati[:idx_values[0]]
                    ti = ti[:idx_values[0]]
            c, o = self.svd_fit(ring, ti, mat=mati, svd_cut=s, apply_correction=apply_correction)
            corr = pandas.concat((corr, c))
            obs = pandas.concat((obs, o))
        return corr, obs
        
        
class LinoptResponseMatrix(ElementResponseMatrix):

    def __init__(self, ring, bpms, quadrupoles, obsnames, obsidx,
                 deltax=1.0e-3, deltay=1.0e-3, **kwargs):
        assert len(obsnames) == len(obsidx), \
            'Linopt RM: observables names and indexes must have the same length' 
        super(LinoptResponseMatrix, self).__init__(okw=kwargs)
        self.set_bpms(ring, bpms, obsnames, obsidx)
        self.set_quadrupoles(ring, quadrupoles, deltax=deltax, deltay=deltay)   
            
    def set_bpms(self, ring, refpts, obsnames, obsidx):
        for o, i in zip(obsnames, obsidx):
            for ii in numpy.atleast_1d(i):
                self.add_observables_refpts(ring, o, refpts=refpts, index=ii)
            
    def set_quadrupoles(self, ring, refpts, deltax=1.0e-3, deltay=1.0e-3):
        self.add_variables_refpts(ring, 'QUAD', deltax, refpts, 'PolynomB', index=1)
        
    def correct(self, ring, target, mat=None, svd_cut=0, apply_correction=False):
        if mat is None:
            mat = self.get_mat()      
        assert len(target) == mat.shape[0], \
            'Linopt RM: target must be of the same length as the matrix'
        corr, obs = self.svd_fit(ring, target, mat=mat, svd_cut=svd_cut, apply_correction=apply_correction)
        return corr, obs                        
