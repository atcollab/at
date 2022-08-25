"""
Classes for response matrix observables and variables
"""

from typing import List
import pandas
import numpy
from at.physics import find_orbit, get_optics


_DATA_NAMES = {'closed_orbit': ['X', 'XP', 'Y', 'YP', 'DP', 'CT'],
               'beta': ['BETX', 'BETY'],
               'alpha': ['ALFX', 'ALFY'],
               'mu': ['MUX', 'MUY'],
               'dispersion': ['DX', 'DPY', 'DY', 'DPY'],
               'W': ['WX', 'WY']}


def sum_polab(ring, refpts, attrname, index):
    k =  [getattr(ring[r], attrname)[index] * 
          getattr(ring[r], 'Length') for r in refpts]                              
    return numpy.sum(k)


def gen_unique_keys(arr):
    vun, vind = numpy.unique(arr, return_inverse=True, axis=0)
    for i in range(len(vun)):
        vr = [ii for ii, vi in enumerate(vind) if vi==i]
        if len(vr) > 1:
            for ii, ik in enumerate(vr):
                arr[ik] = (arr[ik][0], arr[ik][1]+'.'+str(ii))
    return arr


def get_keys(ring, variables, attrkey):
    keys = []
    for v in variables:
        refs = getattr(v, 'refpts', None)
        if refs is None:
            keys.append((v.name, 'GLOBAL'))
        else:
            r = numpy.atleast_1d(refs)[0]               
            keys.append((v.name, getattr(ring[r], attrkey)))
    return pandas.MultiIndex.from_tuples(gen_unique_keys(keys)) 

    
def  load_confs(ring, oconf, vconf):
    vgroups = vconf.index(level=0)
    print(vgroup)
               

      
class Variable(object):

    def __init__(self, name, delta, refpts=None, attname=None, index=None,
                 getsetf=None, args=(), kwargs={}):   
        self.name = name
        self.attname = attname
        self.refpts = refpts
        self.delta = delta
        self.attindex = index  
        self.args = args
        self.kwargs = kwargs
        self.ufun = False
        
        if getsetf is None:   
            assert refpts is not None, \
                'repfts required for element variable'  
            assert attname is not None, \
                'attname required for element variable'                             
            self.setf, self.getf = self._access(index)
        else:
            self.ufun = True
            self.setf, self.getf = getsetf

        super(Variable, self).__init__()     
        
    @staticmethod
    def _access(index):
        """Access to element attributes"""
        if index is None:
            setf = setattr
            getf = getattr
        else:
            def setf(elem, attrname, value):
                getattr(elem, attrname)[index] = value

            def getf(elem, attrname):
                return getattr(elem, attrname)[index]
        return setf, getf       
        
    def set(self, ring, value):
        if self.ufun:
            self.setfun(ring, value, *self.args, **self.kwargs)  
        else:
            for elem in ring.select(self.refpts):
                self.setf(elem, self.attname, value)                      

    def get(self, ring):
        if self.ufun:
            return self.getfun(ring, *self.args, **self.kwargs)
        else:
            values = numpy.array([self.getf(elem, self.attname)
                                  for elem in ring.select(self.refpts)])
            return numpy.average(values)
            
        

class Observable(object):    

    def __init__(self, name, refpts=None, index=None, weight=1, fun=None, args=(), kwargs={}):
        self.refpts = refpts
        self.parname = name
        self.name =name
        self.weight = weight
        self.parindex = index
        self.args = args
        self.kwargs = kwargs
        
        if fun is None:
            assert index is not None, \
                'index required for orbit or linopt observable'
            assert refpts is not None, \
                'repfts required for orbit or linopt observable'
            assert name in _DATA_NAMES.keys(),\
                '{} observable not defined in at.get_optics'.format(name)
            
            self.name = _DATA_NAMES[name][index]
            if self.parname == 'closed_orbit':
                fun = find_orbit
            else:
                fun = get_optics
        self.fun = fun
        self.funname = fun.__name__
        
        super(Observable, self).__init__()
               
    def value(self, ring):
        if self.fun is find_orbit:
            orb0 = self.kwargs.pop('orbit', None)
            _, orbit = self.fun(ring, refpts=self.refpts, orbit=orb0)
            val = orbit[:, self.parindex][0]
        elif self.fun is get_optics:
            get_chrom = (self.parname == 'W')
            twiss_in = self.kwargs.pop('twiss_in', None)
            _, _, ld = self.fun(ring, refpts=self.refpts,
                                get_chrom=get_chrom, twiss_in=twiss_in)
            val = ld[self.parname][self.parindex][0]
        else:
            val = self.fun(ring, *self.args, **self.kwargs)
        assert numpy.isscalar(val), \
            'Observable return value has to be a scalar'
        return val
        
        
class RMElements(List):

    def __init__(self, attrs, attrkey='FamName', attrlist=[], elements=[]):
        self.attrkey = attrkey
        self.attrs = attrs
        self.attrlist = numpy.unique(attrlist+[attrkey])
        self.conf = pandas.DataFrame(columns=numpy.concatenate((attrs,
                                                                attrlist)))
        super(RMElements, self).__init__(elements)   
        
    def add_elements(self, ring, elements):
        self += numpy.atleast_1d(elements) 
        self._update(ring) 
        
    def _update(self, ring):
        keys = get_keys(ring, self, self.attrkey)
        varsd = {} 
        for attr in self.attrs:
            s = [numpy.atleast_1d(getattr(v, attr, None)) for v in self]
            varsd.update({attr: pandas.Series(s, index=keys)})                 
        for attr in self.attrlist:
            s = []
            for v in self:
                ref = numpy.atleast_1d(getattr(v, 'refpts', None))
                if ref[0] is None:
                    s += [ref[0]]
                else:
                    s += [getattr(ring[ref[0]], attr, None)]
            varsd.update({attr: pandas.Series(s, index=keys)})               
        df = pandas.DataFrame(varsd) 
        self.conf = df  
        
        
class RMVariables(RMElements):

    def __init__(self, attrkey='FamName', attrlist=[], variables=[]):
        attrs = ['name', 'refpts', 'attname', 'attindex', 'delta']
        super(RMVariables, self).__init__(attrs, attrkey=attrkey,
                                          attrlist=attrlist,
                                          elements=variables)   
        
    def add_elements_refpts(self, ring, name, delta, refpts, attname, index=None):                          
        attname = numpy.broadcast_to(attname, len(refpts))
        name = numpy.broadcast_to(name, len(refpts))
        delta = numpy.broadcast_to(delta, len(refpts))
        index = numpy.broadcast_to(index, len(refpts))
        self += [Variable(n, d, refpts=r, attname=a, index=i, getsetf=None)
                 for r, a, n, d, i in zip(refpts, attname, name, delta, index)]   
        self._update(ring)        
        
    def get(self, ring, mask=None):
        if mask is None:
            mask = numpy.ones(len(self), dtype=bool)
        return [v.get(ring) for v in numpy.array(self)[mask]]
        
    def set(self, ring, values, mask=None):
        if mask is None:
            mask = numpy.ones(len(self), dtype=bool)
        [v.set(ring, val) for v, val in zip(numpy.array(self)[mask], values)]                  
            
            
class RMObservables(RMElements):

    def __init__(self, attrkey='FamName', attrlist=[], observables=[]):
        self._ref_fun = ['find_orbit', 'get_optics']
        attrs = ['name', 'refpts', 'weight', 'parindex', 'funname']
        super(RMObservables, self).__init__(attrs, attrkey=attrkey,
                                            attrlist=attrlist,
                                            elements=observables)

    def add_elements_refpts(self, ring, name, refpts, index=None, weight=1):
        name = numpy.broadcast_to(name, len(refpts))
        index = numpy.broadcast_to(index, len(refpts))
        weight = numpy.broadcast_to(weight, len(refpts))
        index = numpy.broadcast_to(index, len(refpts))
        self += [Observable(n, r, index=i, weight=w, fun=None)
                 for r, n, i, w in zip(refpts, name, index, weight)]   
        self._update(ring)       

    def values(self, ring, mask=None):
        if mask is None:
            mask = numpy.ones(len(self), dtype=bool)
        allref = [n[0] for n in self.conf.refpts[mask] if n[0] is not None]
        allref = numpy.sort(numpy.unique(allref))
        all_fun = numpy.unique([f[0] for f in self.conf.funname[mask]])
        linopt = 'get_optics' in all_fun
        orbit = 'find_orbit' in all_fun and not linopt
        if linopt:
            get_chrom = numpy.any([n[0]=='W' for n in self.conf.name[mask]])
            _, _, valref = get_optics(ring, refpts=allref, get_chrom=get_chrom)
        if orbit:
            _, valref = find_orbit(ring, refpts=allref)
        vals = []
        for o in numpy.array(self)[mask]:
             if o.refpts is None:
                 vals.append(o.value(ring)) 
             elif linopt:
                 vals.append(valref[allref == o.refpts][o.parname][o.parindex][0])
             elif orbit:
                 vals.append(valref[allref == o.refpts, o.parindex][0])                                   
        return numpy.array(vals)
