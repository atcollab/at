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
import warnings


globring = None
    

class ElementResponseMatrix(object):

    def __init__(self, ring, **kwargs):
        self.ring = ring.copy()
        #  add attributes to ring such that they are
        #  accessible through the global variable
        self.ring.variables = RMVariables(**kwargs)
        self.ring.observables = RMObservables(**kwargs)
        self.fullrm = None
        self.excl_obs = []
        self.excl_var = []
        super(ElementResponseMatrix, self).__init__()
        
    def add_variables(self, variables):
        self.ring.variables.add_elements(self.ring, variables)
        
    def add_observables(self, observables):
        self.ring.observables.add_elements(self.ring, observables)
        
    def add_variables_refpts(self, name, delta, refpts, attname,
                             index=None, sum_zero=False):
        self.ring.variables.add_elements_refpts(self.ring, name, 
                                                delta, refpts,
                                                attname, index=index)
        if sum_zero:
            assert 'Polynom' in attname,\
               'sum_zero available only for PolynomA/B attribute'
            assert index is not None,\
               'index required for sum_zero'
            sum_zero = Observable(name+'_SUM', fun=sum_polab,
                                  args=(refpts, attname, index))
            self.add_observables(sum_zero)
                                                
    def add_observables_refpts(self, name, refpts, index=None, weight=1):
        self.ring.observables.add_elements_refpts(self.ring, name,
                                                  refpts, index=index,
                                                  weight=weight)
                                  
    @staticmethod
    def _resp_one(ring, variable, key):
        if ring is None:
            ring = globring
        v0 = variable.get(ring)
        variable.set(ring, v0+variable.delta)
        op = ring.observables.values(ring)
        variable.set(ring, v0-variable.delta)  
        om = ring.observables.values(ring)
        do = {key: (oom-oop)/2/variable.delta for key, oom, oop 
              in zip(ring.observables.conf.index, om, op)}
        variable.set(ring, v0)
        return {key: pandas.Series(do)}         
        
    def compute_fullrm(self, use_mp=False, pool_size=None, start_method=None):
        rv = {}
        if use_mp:
            ctx = multiprocessing.get_context(start_method)
            if pool_size == None:
                pool_size = min(len(self.ring.variables),
                                multiprocessing.cpu_count())
            if ctx.get_start_method() == 'fork':
                global globring
                globring = self.ring
                args = [(None, var, key) for var, key in zip(self.ring.variables,
                                                             self.ring.variables.conf.index)]
                with ctx.Pool(pool_size) as pool:
                    results = pool.starmap(partial(self._resp_one), args)
                globring = None
            else:
                args = [(self.ring, var, key) for var in zip(self.ring.variables,
                                                             self.ring.variables.conf.index)]
                with ctx.Pool(pool_size) as pool:
                    results = pool.starmap(partial(self._resp_one), args)
            for r in results:
                rv.update(r)             
        else: 
            for var, key in zip(self.ring.variables, self.ring.variables.conf.index):  
                rv.update(self._resp_one(self.ring, var, key))
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
            store['vconf'] = self.ring.variables.conf
            store['oconf'] = self.ring.observables.conf
        store.close()
        
    def load_fullrm(self, filename):
        store = pandas.HDFStore(filename)
        self.fullrm = store['fullrm']
        vconf = store['vconf']
        oconf = store['oconf']
        store.close()
        load_confs(self, oconf, vconf)
        
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
            self.excl_var = remove_element(vardictlist, self.ring.variables.conf)
        if obsdictlist is not None:
            self.excl_obs = remove_element(obsdictlist, self.ring.observables.conf)
       
    def get_mat(self):
        assert self.fullrm is not None,\
           ' Empty response matrix: please run compute_fullrm() first'                               
        return self.fullrm.loc[~self.fullrm.index.isin(self.excl_obs),
                               ~self.fullrm.columns.isin(self.excl_var)]                           
        
    def get_vals(self, mat=None):
        if mat is None:
            mat = self.get_mat()
        mask = self.ring.observables.conf.index.isin(mat.index)  
        obs = numpy.array(self.ring.observables.values(self.ring, mask))
        mask = self.ring.variables.conf.index.isin(mat.columns)
        var = numpy.array(self.ring.variables.get(self.ring, mask))
        return obs, var
        
    def svd(self, mat=None):
        if mat is None:
            mat = self.get_mat()
        return numpy.linalg.svd(mat, full_matrices=False)                                                
                
    def svd_fit(self, target, mat=None, svd_cut=0, apply_correction=False):
        if mat is None:
            mat = self.get_mat()
        val, var = self.get_vals(mat=mat)
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
            mask = self.ring.variables.conf.index.isin(mat.columns) 
            self.ring.variables.set(self.ring, var+dk, mask)
        valn, varn = self.get_vals(mat=mat) 
        corr.update({'New': pandas.Series(varn, mat.columns)})
        obs.update({'New': pandas.Series(valn, mat.index)})
        return pandas.DataFrame(corr), pandas.DataFrame(obs)

        
class OrbitResponseMatrix(ElementResponseMatrix):

    _MAT_INDEX = {'X': ['X', 'HST_SUM'],
                  'Y': ['Y', 'VST_SUM']}
    _MAT_COLS = {'X': ['HST', 'RF'], 'Y': ['VST']}         

    def __init__(self, ring, bpms, steerers, cavities=None, sum_zero=True):
        super(OrbitResponseMatrix, self).__init__(ring)
        self.set_bpms(bpms)
        self.set_steerers(steerers, sum_zero=sum_zero)
        if cavities is None:
            ring = ring.radiation_off(copy=True)
        else:
            ring = ring.radiation_off(cavity_pass='RFCavityPass', copy=True)
            self.set_cavities(cavities)   
        
    def set_bpms(self, refpts):
        self.add_observables_refpts('closed_orbit', refpts=refpts, index=0)
        self.add_observables_refpts('closed_orbit', refpts=refpts, index=2)
        
    def set_steerers(self, refpts, sum_zero=True):
        self.add_variables_refpts('HST', 1.0e-4, refpts,'PolynomB', index=0, sum_zero=sum_zero)
        self.add_variables_refpts('VST', 1.0e-4, refpts,'PolynomA', index=0, sum_zero=sum_zero)
        
    def set_cavities(self, refpts):
        self.add_variables_refpts('RF', 10, [refpts],'Frequency')
        
    def set_excluded_bpms(self, bpmdictlist):
        self.set_exluded(obsdistlist=bpmdictlist)
        
    def set_excluded_steerer(self, steerersdictlist):
        self.set_exluded(vardistlist=steerersdictlist)
        
    def correct_orbit(self, target=0, plane=['X', 'Y'], svd_cut=0, apply_correction=False):
        mat = self.get_mat()
        svd_cut = numpy.broadcast_to(svd_cut, len(plane))
        target = numpy.broadcast_to(target, len(plane))
        corr = pandas.DataFrame()
        obs = pandas.DataFrame()
        for t, p, s in zip(target, plane, svd_cut):
            assert p=='X' or p=='Y',\
               'Orbit plane has to be X or Y'
            mati = mat.loc[mat.index.isin(self._MAT_INDEX[p], level=0),
                           mat.columns.isin(self._MAT_COLS[p], level=0)]
            c, o = self.svd_fit(t, mat=mati, svd_cut=s, apply_correction=apply_correction)
            corr = pandas.concat((corr, c))
            obs = pandas.concat((obs, o))
        return corr, obs                             
