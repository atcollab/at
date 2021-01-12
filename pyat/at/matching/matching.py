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
        return '\n{:>12s}{:>13s}{:>16s}{:>16s}\n'.format(
            'Name', 'Initial', 'Final ', 'Variation')

    def status(self, ring, vini=np.NaN):
        vnow = self.get(ring)
        return '{:>12s}{: 16e}{: 16e}{: 16e}'.format(
            self.name, vini, vnow, (vnow - vini) / vini)


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
          cnstrs.add(circ_fun, 850.0)
          cnstrs.add(mcf_fun, 1.0e-4, weight=0.1)
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

    def add(self, fun, target, name=None, weight=1.0, bounds=(0, 0)):
        """Add a target to the Constraints container

        PARAMETERS
            fun           evaluation function. Called as:
                          value = fun(ring, *args, **kwargs)
                            value is the constrained parameter value
                            value may be a scalar or an array.
                            the positional and keyword parameters come from
                            the Constraints initialisation
            target        desired value.

        KEYWORDS
            name=None     name of the constraint. If None, name is generated
                          from the name of the evaluation function
            weight=1.0    weight factor: the residual is (value-target)/weight
            bounds=(0,0)  lower and upper bounds. The parameter is constrained
                          in the interval [target-low_bound target+up_bound]

        The "target", "weight" and "bounds" input must be broadcastable to the
        shape of "value".
        """
        if name is None:                # Generate the constraint name
            name = fun.__name__
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
            weight = np.broadcast_to(weight, value.shape)
            diff = value - target
            lb = np.ravel(diff - lbound)
            ub = np.ravel(diff - ubound)
            lb[lb >= 0] = 0
            ub[ub <= 0] = 0
            return np.maximum(abs(lb), abs(ub)) / np.ravel(weight)

        return np.concatenate([res(v, t, w, lb, ib) for v, t, w, lb, ib
                               in zip(self.values(ring), self.target,
                                      self.weight, self.lbound,
                                      self.ubound)], axis=None)

    @staticmethod
    def header():
        """Header for the display of constraint values and residuals"""
        return '\n{:>12s}{:>13s}{:>16s}{:>16s}{:>16s}\n'.format(
            'Name', 'Initial', 'Final ', 'Target', 'Residual')

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
                strs.append('{:>12s}{: 16e}{: 16e}{: 16e}{: 16e}'.format(
                    name, vi, vn, vt, vn - vt))
        return '\n'.join(strs)


class _RefConstraints(Constraints):
    """Base class for position-related constraints: handle the refpoints
    of each target"""
    def __init__(self, ring, *args, **kwargs):
        self.nelems = len(ring)
        self.refs = []
        self.refpts = bool_refpts([], self.nelems)
        super(_RefConstraints, self).__init__(*args, **kwargs)

    def add(self, fun, target, refpts=None, **kwargs):
        ref = bool_refpts(refpts, self.nelems)
        # Store the new refpoint
        self.refs.append(ref)
        # Update the union of all refpoints
        self.refpts = np.stack((self.refpts, ref), axis=0).any(axis=0)
        super(_RefConstraints, self).add(fun, target, **kwargs)

    def values(self, ring):
        # Single optics computation
        vloc, vglob = self.compute(ring, *self.args, **self.kwargs)
        # Evaluate all constraints
        return [fun(loc, *vglob) for fun, loc in zip(self.fun, vloc)]

    def compute(self, ring, *args, **kwargs):
        """Dummy computation. Compute must return:
        - an iterator over local data for each target
        - a tuple of global data"""
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
          cnstrs.add('tune', 0.44, order=0, weight=0.01)

          # Add a chromaticity constraint (both planes)
          cnstrs.add('chrom', [0.0 0.0])

          # define a constraint of phase advance between 2 points
          def mu_diff(lindata, tune, chrom):
              delta_mu = (lindata[1].mu[0] - lindata[0].mu[0])/(2*np.pi)
              return delta_mu % 1.0

          # Add a phase advance constraint, giving the desired locations
          cnstrs.add(mu_diff, 0.5, refpts=[sf0 sf1])
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

    def add(self, param, target, refpts=None, order=None, name=None, **kwargs):
        """Add a target to the LinoptConstraints container

        PARAMETERS
            param         2 possibilities:
                          - parameter name: see at.linopt for the name of
                            available parameters. In addition to local optical
                            parameters, 'tune' and 'chrom' are allowed.
                          - user-supplied parameter evaluation function:
                                value = param(lindata, tune, chrom)
                            lindata contains the optics parameters at all the
                              specified refpoints
                            value is the constrained parameter value
                              value may be a scalar or an array.
            target        desired value.

        KEYWORDS
            refpts=None   location of the constraint. Several locations may be
                          given to apply the same constraint at several points.
            order=None    index in the parameter array. If None, the full array
                          is used.
            name=None     name of the constraint. If None, name is generated
                          from param and order.
            weight=1.0    weight factor: the residual is (value-target)/weight.
            bounds=(0,0)  lower and upper bounds. The parameter is constrained
                          in the interval [target-low_bound target+up_bound]

        The target, weight and bounds values must be broadcastable to the shape
        of value.
        """
        # noinspection PyUnusedLocal
        def tunefun(refdata, tune, chrom):
            return getv(tune)

        # noinspection PyUnusedLocal
        def chromfun(refdata, tune, chrom):
            return getv(chrom)

        # noinspection PyUnusedLocal
        def attrfun(refdata, tune, chrom):
            return getf(refdata, param)

        getf, getv = _dataaccess(order)

        if name is None:                # Generate the constraint name
            if callable(param):
                name = param.__name__
            elif order is None:
                name = param
            else:
                name = '{}_{}'.format(param, order)

        if callable(param):
            fun = param
            # self.refpts[:] = True     # necessary not to miss 2*pi jumps
            self.get_chrom = True       # fun may use dispersion or chroma
        elif param == 'tune':
            fun = tunefun
            refpts = []
        elif param == 'chrom':
            fun = chromfun
            refpts = []
            self.get_chrom = True       # slower but necessary
        else:
            fun = attrfun
            if param == 'dispersion':
                self.get_chrom = True   # slower but necessary
            # elif param == 'mu':
            #     self.refpts[:] = True # necessary not to miss 2*pi jumps

        super(LinoptConstraints, self).add(fun, target, refpts, name=name,
                                           **kwargs)

    def compute(self, ring, *args, **kwargs):
        """Optics computation before evaluation of all constraints"""
        ld0, tune, chrom, ld = linopt(ring, refpts=self.refpts,
                                      get_chrom=self.get_chrom, **kwargs)
        return (ld[ref[self.refpts]] for ref in self.refs), (tune, chrom)


def match(ring, variables, constraints, verbose=2, max_nfev=1000,
                  diff_step=1.0e-10):
    """Perform matching of constraints by varying variables

    PARAMETERS
        ring                ring lattice or transfer line
        variables           sequence of Variable objects
        constraints         sequance of Constraints objects

    KEYWORDS
        verbose=2           See scipy.optimize.least_squares
        max_nfev=1000           "
        diff_step=1.0e-10       "
    """
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
    least_squares(fun, init, bounds=bounds, verbose=verbose, max_nfev=max_nfev,
                  method=method, diff_step=diff_step)

    print(Constraints.header())
    for cst, ini in zip(constraints, initvals):
        print(cst.status(ring, initial=ini))

    print(Variable.header())
    for var, vini in zip(variables, init):
        print(var.status(ring, vini=vini))
    print()
