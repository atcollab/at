import numpy as np
from scipy.optimize import least_squares
from itertools import repeat
from at.lattice import refpts_iterator, bool_refpts
from at.physics import linopt


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


def _dataaccess(index):
    """Access to optics parameters"""
    if index is None:
        getf = getattr

        def getv(x):
            return x
    else:
        def getf(elem, attrname):
            return getattr(elem, attrname)[:, index]

        def getv(x):
            return x[index]
    return getf, getv


class Variable(object):
    """A Variable is a scalar value acting on a lattice through the
    user-defined functions setfun and getfun
    """
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

    @staticmethod
    def header():
        return '{:s}\t{:s}\t\t{:s}\t\t{:s}'.format('Name', 'Initial', 'Final',
                                                   'Variation')

    def status(self, ring, vini=np.NaN):
        vnow = self.get(ring)
        return '{:s}\t{:e}\t{:e}\t{:e}'.format(self.name, vini, vnow,
                                               (vnow - vini) / vini)


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
    """Container for generic constraints:
      - a constraint is defined by a user-defined evaluation function.
      - constraints are added to the container with the Constraints.add method

      Example:
          # define an evaluation function for the ring circumference:
          def circ_fun(ring):
              return ring.get_s_pos(len(ring) + 1)

          # define an evaluation function for the momentum compaction factor:
          def mcf_fun(ring):
              return ring.get_mcf()

          # Construct the container:
          cnstrs = Constraints()

          # Add the two constraints:
          cnstrs.add('circ', circ_fun, 850.0)
          cnstrs.add('mcf', mcf_fun, 1.0e-4, weight=0.1)
    """
    def __init__(self, *args, **kwargs):
        """Constraints(*args, **kwargs)
        build a generic constraints container.

        The positional and keyword parameters are provided to all the evaluation
        functions.
        """
        self.name = []
        self.fun = []
        self.target = []
        self.weight = []
        self.lbound = []
        self.ubound = []
        self.args = args
        self.kwargs = kwargs

    def add(self, name, fun, target, weight=1.0, bounds=(0, 0)):
        """Add a target to the Constraints container

        PARAMETERS
            name          name of the constraint
            fun           evaluation function. Called as:
                          value = fun(ring, *args, **kwargs)
                            value is the constrained parameter value
                            value may be a scalar or an array.
                            the positional and keyword parameters come from
                            the Constraints initialisation
            target        desired value.

        KEYWORDS
            weight=1.0    weight factor: the residual is (value-target)/weight
            bounds=(0,0)  lower and upper bounds with respect to target

        The target, weight and bounds values must be broadcastable to the shape
        of value.
        """
        self.name.append(name)
        self.fun.append(fun)
        self.target.append(target)
        self.weight.append(weight)
        self.lbound.append(bounds[0])
        self.ubound.append(bounds[1])

    def values(self, ring):
        """Return the list of actual parameter values"""
        return [fun(ring, *self.args, **self.kwargs) for fun in self.fun]

    def evaluate(self, ring):
        """Return a flattened array of weighted residuals"""

        def res(value, target, weight, lbound, ubound):
            diff = value - target
            lb = diff - lbound
            ub = diff - ubound
            lb[lb >= 0] = 0
            ub[ub <= 0] = 0
            return np.maximum(abs(lb), abs(ub)) / weight

        return np.concatenate([res(v, t, w, lb, ib) for v, t, w, lb, ib
                               in zip(self.values(ring), self.target,
                                      self.weight, self.lbound,
                                      self.ubound)], axis=None)

    @staticmethod
    def header():
        """Header for the display of constraint residuals"""
        return '{:s}\t{:s}\t\t{:s}\t\t{:s}\t\t{:s}'.format('Name', 'Initial',
                                                           'Final', 'Target',
                                                           'Residual')

    def status(self, ring, initial=None):
        """Return a string giving the actual state of constraints"""
        if initial is None:
            initial = repeat(np.NaN)
        strs = []
        for name, ini, now, target in zip(self.name, initial,
                                          self.values(ring), self.target):
            for vi, vn, vt in zip(np.broadcast_to(ini, now.shape).flat,
                                  now.flat,
                                  np.broadcast_to(target, now.shape).flat):
                strs.append('{:s}\t{:e}\t{:e}\t{:e}\t{:e}'.
                            format(name, vi, vn, vt, vn - vt))
        return '\n'.join(strs)


class _RefConstraints(Constraints):
    """Base class for position-related constraints"""
    def __init__(self, ring, *args, **kwargs):
        self.nelems = len(ring)
        self.refs = []
        self.refpts = bool_refpts([], self.nelems)
        super(_RefConstraints, self).__init__(*args, **kwargs)

    # noinspection PyMethodOverriding
    def add(self, name, refpts, fun, target, weight=1.0, bounds=(0, 0)):
        ref = bool_refpts(refpts, self.nelems)
        self.refs.append(ref)
        self.refpts = np.stack((self.refpts, ref), axis=0).any(axis=0)
        super(_RefConstraints, self).add(name, fun, target, weight, bounds)

    def values(self, ring):
        vloc, vglob = self.compute(ring, *self.args, **self.kwargs)
        return [fun(loc, *vglob) for fun, loc in zip(self.fun, vloc)]

    def compute(self, ring, *args, **kwargs):
        return repeat(None), ()


class LinoptConstraints(_RefConstraints):
    """Container for linear optics constraints:
      - a constraint can be set on any result of at.linopt
      - constraints are added to the container with the LinoptConstraints.add
        method.

      at.linopt is called once before the evaluation of all constraints

      Example:
          cnstrs = LinoptConstraints(ring, dp=0.01, coupled=False)

          # Add a beta H (beta[0]) constraint at location ref_inj
          cnstrs.add('beta_x_inj', 'beta', 18.0, refpts=ref_inj, order=0)

          # Add a tune constraint
          cnstrs.add('tune_x', 'tune', 0.44, order=0, weight=0.01)

          # Add a chromaticity constraint 9both planes
          cnstrs.add('chromaticity', 'chrom', [0.0 0.0])
    """
    def __init__(self, ring, *args, **kwargs):
        """Build a LinoptConstraints container

        KEYWORDS
        dp=0.0          momentum deviation.
        XYStep=1.0e-8   transverse step for numerical computation
        DPStep=1.0E-6   momentum deviation used for computation of
                        chromaticities and dispersion
        coupled=True    if False, simplify the calculations by assuming
                        no H/V coupling
        twiss_in=None   Initial twiss parameters for transfer line optics.
                        "lindata" stucture, where only the beta and alpha are
                        required and used.
        orbit=None      Initial trajectory for transfer line
                        ((6,) array)
        """
        self.get_chrom = False
        super(LinoptConstraints, self).__init__(ring, *args, **kwargs)

    def add(self, name, attrname, target, refpts=None, order=None, weight=1.0,
            bounds=(0, 0)):
        """Add a target to the LinoptConstraints container

        PARAMETERS
            name          name of the constraint
            attrname      parameter name: see at.linopt for the name of
                          available parameters. In addition to local optical
                          parameters, 'tune' and 'chrom' are allowed.
            target        desired value.

        KEYWORDS
            refpts=None   location of the constraint
            order=None    index in the parameter array. If None, the full array
                          is used.
            weight=1.0    weight factor: the residual is (value-target)/weight
            bounds=(0,0)  lower and upper bounds with respect to target

        The target, weight and bounds values must be broadcastable to the shape
        of value.
        """
        getf, getv = _dataaccess(order)

        # noinspection PyUnusedLocal
        def tunefun(refdata, tune, chrom):
            return getv(tune)

        # noinspection PyUnusedLocal
        def chromfun(refdata, tune, chrom):
            return getv(chrom)

        # noinspection PyUnusedLocal
        def attrfun(refdata, tune, chrom):
            return getf(refdata, attrname)

        if attrname == 'tune':
            fun = tunefun
            refpts = []
        elif attrname == 'chrom':
            fun = chromfun
            refpts = []
            self.get_chrom = True
        else:
            fun = attrfun
            if attrname == 'dispersion':
                self.get_chrom = True

        super(LinoptConstraints, self).add(name, refpts, fun, target, weight,
                                           bounds)

    def compute(self, ring, *args, **kwargs):
        """Optics computation before evaluation on all constraints"""
        ld0, tune, chrom, ld = linopt(ring, refpts=self.refpts,
                                      get_chrom=self.get_chrom, **kwargs)
        return (ld[ref[self.refpts]] for ref in self.refs), (tune, chrom)


class oldConstraints(object):
    # noinspection PyUnusedLocal
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


def match(ring, variables, constraints):
    def fun(vals):
        for value, variable in zip(vals, variables):
            variable.set(ring, value)
        return np.concatenate([cons.evaluate(ring) for cons in constraints],
                              axis=None)

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

    initvals = [cst.values(ring) for cst in constraints]
    least_squares(fun, init, bounds=bounds, verbose=2, max_nfev=1000,
                  method=method, diff_step=1.0e-10)

    print()
    print(Constraints.header())
    print()
    for cst, ini in zip(constraints, initvals):
        print(cst.status(ring, initial=ini))
    print()
    print(Variable.header())
    for var, vini in zip(variables, init):
        print(var.status(ring, vini=vini))
    print()
