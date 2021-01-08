import numpy as np
from scipy.optimize import least_squares
from at.lattice import refpts_iterator
from at.physics import linopt


def _access(order):
    if order is None:
        setf = setattr
        getf = getattr
    else:
        def setf(elem, attrname, value):
            getattr(elem, attrname)[order] = value

        def getf(elem, attrname):
            return getattr(elem, attrname)[order]
    return setf, getf


class Variable(object):
    """A Variable is a scalar value acting on a lattice through the
    user-defined functions setfun and getfun"""

    def __init__(self, setfun, getfun, name='', bounds=(-np.inf, np.inf),
                 args=()):
        self.setfun = setfun
        self.getfun = getfun
        self.name = name
        self.bounds = bounds
        self.args = args
        super(Variable, self).__init__()

    def set(self, ring, value):
        self.setfun(ring, value, *self.args)

    def get(self, ring):
        return self.getfun(ring, *self.args)


class ElementVariable(Variable):
    """An ElementVariable is:
    - a scalar attribute or
    - an element of an array attribute
    of one or several elements of a lattice"""

    def __init__(self, refpts, attname, name='', order=None,
                 bounds=(-np.inf, np.inf)):
        setf, getf = _access(order)

        def setfun(ring, value):
            for elem in refpts_iterator(ring, refpts):
                setf(elem, attname, value)

        def getfun(ring):
            values = np.array([getf(elem, attname) for elem in
                               refpts_iterator(ring, refpts)])
            return np.average(values)

        super(ElementVariable, self).__init__(setfun, getfun, name, bounds)


class Constraints(object):
    def __init__(self, ring, fun, constraints, args=()):
        self.fun = fun
        self.constraints = constraints
        self.args = args
        self.name = np.array([ci['name'] for ci in self.constraints])
        self.target = np.array([ci['target'] for ci in self.constraints])
        self.bounds = np.array([ci['bounds'] for ci in self.constraints]).T
        self.weights = np.array([ci['weight'] for ci in self.constraints])

    @staticmethod
    def build(varname, target, bounds=(0, 0), weight=1):
        return {'name': varname, 'target': target,
                'bounds': bounds, 'weight': weight}

    def get_vals(self, ring):
        return self.fun(ring, *self.args)

    def evaluate(self, ring):
        diff = self.get_vals(ring) - self.target
        lb = diff - self.bounds[0, :]
        ub = diff - self.bounds[1, :]
        lb[lb >= 0] = 0
        ub[ub <= 0] = 0
        diff = np.maximum(abs(lb), abs(ub)) * self.weights
        return diff


class RingConstraints(object):
    def __init__(self, ring, refpts, constraints, args=()):
        self.refpts = refpts
        self.constraints = constraints
        self.args = args
        self.name, self.target, self.bounds = self.get_celems(ring)

    @staticmethod
    def build(varname, elemt, target, bounds=(0, 0)):
        return {'name': varname, 'elemt': elemt, 'target': target,
                'bounds': bounds}

    def get_celems(self, ring):
        varnames = []
        target = []
        bounds = []
        for r, c in zip(self.refpts, self.constraints):
            ename = ring[r].FamName
            for ci in c:
                varnames.append(
                    ename + '_' + ci['name'] + '_' + str(ci['elemt']))
                target.append(ci['target'])
                bounds.append(ci['bounds'])
        return np.array(varnames), np.array(target), np.array(bounds).T

    def get_vals(self, ring):
        l0, tune, chrom, ld = linopt(ring, refpts=self.refpts, *self.args)
        return np.array([ld[i][ci['name']][ci['elemt']] for i, c in
                         enumerate(self.constraints) for ci in c])

    def evaluate(self, ring):
        diff = self.get_vals(ring) - self.target
        lb = diff - self.bounds[0, :]
        ub = diff - self.bounds[1, :]
        lb[lb >= 0] = 0
        ub[ub <= 0] = 0
        diff = np.maximum(abs(lb), abs(ub))
        return diff


class MultiConstraints(object):
    def __init__(self, constraints):
        self.constraints = constraints
        self.name = self.constraints[0].name
        self.target = self.constraints[0].target
        self.bounds = self.constraints[0].bounds
        for c in self.constraints[1:]:
            self.name = np.concatenate((self.name, c.name))
            self.target = np.concatenate((self.target, c.target))
            self.bounds = np.concatenate((self.bounds, c.bounds))

    def get_vals(self, ring):
        vals = []
        for c in self.constraints:
            vals = np.concatenate((vals, c.get_vals(ring)))
        return np.array(vals)

    def evaluate(self, ring):
        diff = []
        for c in self.constraints:
            diff = np.concatenate((diff, c.evaluate(ring)))
        return np.array(diff)


def match(ring, variables, constraints):
    def fun(vals, constraints, variables):
        for val, variable in zip(vals, variables):
            variable.set(ring, val)
        return constraints.evaluate(ring)

    init = []
    bounds = []
    for var in variables:
        init.append(var.get(ring))
        bounds.append(var.bounds)
    bounds = np.squeeze(bounds).T

    if np.all(bounds == np.inf) and np.size(constraints.target) >= np.size(
            variables):
        method = 'lm'
    else:
        method = 'trf'
    print(' ')
    print('Using method', method)
    print(' ')

    initvals = constraints.get_vals(ring)
    least_squares(fun, init, bounds=bounds, args=(constraints, variables),
                  verbose=2, max_nfev=1000, method=method, diff_step=1.0e-10)
    finalvals = constraints.get_vals(ring)
    print(' ')
    print(
        '{:s}\t{:s}\t\t{:s}\t\t{:s}\t\t{:s}'.format('Name', 'Initial', 'Final',
                                                    'Target', 'Residual'))
    for iv, fv, name, target in zip(initvals, finalvals, constraints.name,
                                    constraints.target):
        print('{:s}\t{:e}\t{:e}\t{:e}\t{:e}'.format(name, iv, fv, target,
                                                    (fv - target)))
    print(' ')
    print('{:s}\t{:s}\t\t{:s}\t\t{:s}'.format('Name', 'Initial', 'Final',
                                              'Variation'))
    for i, var in enumerate(variables):
        v0 = init[i]
        v = var.get(ring)
        print('{:s}\t{:e}\t{:e}\t{:e}'.format(var.name, v0, v, (v - v0) / v0))
    print(' ')
