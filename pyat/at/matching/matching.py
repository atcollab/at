from scipy.optimize import minimize, least_squares 
import at
import numpy as np


class FunVariable(object):
    def __init__(self,ring,vfun,gfun,varname='',bounds=(-np.inf,np.inf),args=()):
        self.vfun = vfun
        self.gfun = gfun
        self.varname = varname
        self.bounds=bounds
        self.args = args
        self.val0 = np.average(self.get_val(ring))
        super(FunVariable, self).__init__()

    def vary(self,ring,val):     
        self.vfun(val,ring,*self.args)

    def get_val(self,ring):
        return self.gfun(ring,*self.args)


class RingVariable(object):
    def __init__(self,ring,refpts,attname,varname='',order=None,bounds=(-np.inf,np.inf)):
        self.refpts = refpts
        self.attname=attname
        self.varname = varname
        self.order=order
        self.bounds=bounds
        self.val0 = np.average(self.get_val(ring))
        super(RingVariable, self).__init__()

    def vary(self,ring,val):      
        at.set_value_refpts(ring, self.refpts, self.attname, val, 
                            order=self.order, increment=False)

    def get_val(self,ring):
        val = at.get_value_refpts(ring, self.refpts, self.attname, 
                                  order=self.order)
        return val


class Constraints(object):
    def __init__(self,ring, fun, constraints, args=()):
        self.fun = fun
        self.constraints = constraints
        self.args = args
        self.varname = np.array([ci['varname'] for ci in self.constraints])
        self.target =  np.array([ci['target'] for ci in self.constraints])
        self.bounds =  np.array([ci['bounds'] for ci in self.constraints]).T
        self.weights =  np.array([ci['weight'] for ci in self.constraints])

    @staticmethod
    def build(varname,target,bounds=(0,0),weight=1):
        return {'varname': varname,'target':target,
                'bounds':bounds,'weight':weight}

    def get_vals(self,ring):
        return self.fun(ring, *self.args)

    def evaluate(self,ring):
        diff = self.get_vals(ring) - self.target
        lb = diff - self.bounds[0,:] 
        ub = diff - self.bounds[1,:]
        lb[lb>=0]=0
        ub[ub<=0]=0
        diff = np.maximum(abs(lb),abs(ub))*self.weights
        return diff


class RingConstraints(object):
    def __init__(self,ring, refpts, constraints, args=()):
        self.refpts = refpts
        self.constraints = constraints
        self.args = args
        self.varname, self.target, self.bounds = self.get_celems(ring)

    @staticmethod
    def build(varname,elemt,target,bounds=(0,0)):
        return {'varname': varname,'elemt':elemt,'target':target,'bounds':bounds}

    def get_celems(self,ring):
        varnames = []
        target = []
        bounds = []
        for r,c in zip(self.refpts,self.constraints):
            ename = ring[r].FamName
            for ci in c:
                varnames.append(ename+'_'+ci['varname']+'_'+str(ci['elemt']))
                target.append(ci['target'])
                bounds.append(ci['bounds'])
        return np.array(varnames),np.array(target),np.array(bounds).T
             
    def get_vals(self,ring):
        l0, tune, chrom, ld = at.linopt(ring,refpts=self.refpts, *self.args)
        return np.array([ld[i][ci['varname']][ci['elemt']] for i,c in enumerate(self.constraints) for ci in c])

    def evaluate(self,ring):
        diff = self.get_vals(ring) - self.target
        lb = diff - self.bounds[0,:] 
        ub = diff - self.bounds[1,:]
        lb[lb>=0]=0
        ub[ub<=0]=0
        diff = np.maximum(abs(lb),abs(ub))
        return diff


class MultiConstraints(object):
    def __init__(self,constraints):
        self.constraints = constraints 
        self.varname = self.constraints[0].varname
        self.target = self.constraints[0].target
        self.bounds = self.constraints[0].bounds
        for c in self.constraints[1:]:
           self.varname = np.concatenate((self.varname,c.varname))
           self.target = np.concatenate((self.target,c.target))
           self.bounds = np.concatenate((self.bounds,c.bounds))
        
    def get_vals(self,ring):
        vals = []
        for c in self.constraints:
            vals = np.concatenate((vals,c.get_vals(ring))) 
        return np.array(vals)     

    def evaluate(self,ring):
        diff = []
        for c in self.constraints:
            diff = np.concatenate((diff,c.evaluate(ring)))
        return np.array(diff)   
            
                            
def match(ring,variables,constraints):

    def fun(vals,constraints,variables):
        for val,variable in zip(vals,variables):
            variable.vary(ring,val)
        return constraints.evaluate(ring)

    init = []
    bounds = []
    for var in variables:
       init.append(var.val0)
       bounds.append(var.bounds)
    bounds = np.squeeze(bounds).T

    if np.all(bounds==np.inf) and np.size(constraints.target)>=np.size(variables):
        method = 'lm'
    else:
        method = 'trf'
    print(' ')
    print('Using method',method)    
    print(' ') 

    initvals  = constraints.get_vals(ring)
    least_squares(fun, init, bounds=bounds,args=(constraints,variables),verbose=2,max_nfev= 1000,method=method,diff_step=1.0e-10)
    finalvals = constraints.get_vals(ring)
    print(' ')
    print('{:s}\t{:s}\t\t{:s}\t\t{:s}\t\t{:s}'.format('Name','Initial','Final','Target','Residual'))
    for iv,fv,name,target in zip(initvals,finalvals,constraints.varname,constraints.target):
        print('{:s}\t{:e}\t{:e}\t{:e}\t{:e}'.format(name,iv,fv,target,(fv-target)))      
    print(' ')
    print('{:s}\t{:s}\t\t{:s}\t\t{:s}'.format('Name','Initial','Final','Variation'))
    for var in variables:
        val = np.average(var.get_val(ring))
        print('{:s}\t{:e}\t{:e}\t{:e}'.format(var.varname,var.val0,val,(val-var.val0)/var.val0))
    print(' ')
