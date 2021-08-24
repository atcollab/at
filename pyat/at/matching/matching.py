from itertools import chain
import numpy as np
from scipy.optimize import least_squares
from itertools import repeat
from at.lattice import refpts_iterator, bool_refpts, uint32_refpts
from at.physics import get_optics, ohmi_envelope, find_orbit


class Variable(object):
    """A Variable is a scalar value acting on a lattice through the
    user-defined functions setfun and getfun
    """
    def __init__(self, setfun, getfun, name='', bounds=(-np.inf, np.inf),
                 *args, **kwargs):
        self.setfun = setfun
        self.getfun = getfun
        self.name = name
        self.bounds = bounds
        self.args = args
        self.kwargs = kwargs
        super(Variable, self).__init__()

    def set(self, ring, value):
        self.setfun(ring, value, *self.args, **self.kwargs)

    def get(self, ring):
        return self.getfun(ring, *self.args, **self.kwargs)

    @staticmethod
    def header():
        return '\n{:>12s}{:>13s}{:>16s}{:>16s}\n'.format(
            'Name', 'Initial', 'Final ', 'Variation')

    def status(self, ring, vini=np.NaN):
        vnow = self.get(ring)
        return '{:>12s}{: 16e}{: 16e}{: 16e}'.format(
            self.name, vini, vnow, (vnow - vini))


class ElementVariable(Variable):
    """An ElementVariable is:
    - a scalar attribute or
    - an element of an array attribute
    of one or several elements of a lattice"""

    def __init__(self, refpts, attname, index=None, **kwargs):
        setf, getf = self._access(index)

        def setfun(ring, value):
            for elem in refpts_iterator(ring, refpts):
                setf(elem, attname, value)

        def getfun(ring):
            values = np.array([getf(elem, attname) for elem in
                               refpts_iterator(ring, refpts)])
            return np.average(values)

        super(ElementVariable, self).__init__(setfun, getfun, **kwargs)
        self.refpts = refpts

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
        self.rad = kwargs.pop('rad', False)
        self.args = args
        self.kwargs = kwargs

    def add(self, fun, target, name=None, weight=1.0, bounds=(0.0, 0.0)):
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
        target, lbound, ubound = np.atleast_1d(target, *bounds)
        self.target.append(target)
        self.weight.append(weight)
        self.lbound.append(lbound)
        self.ubound.append(ubound)

    def values(self, ring):
        """Return the list of actual parameter values"""
        return [fun(ring, *self.args, **self.kwargs) for fun in self.fun]

    def evaluate(self, ring):
        """Return a flattened array of weighted residuals"""

        def res(value, target, weight, lbound, ubound):
            diff = value - target   # broadcast here => at least 1 dimension
            lb = diff + lbound
            ub = diff - ubound
            lb[lb >= 0] = 0         # needs at least 1 dimension
            ub[ub <= 0] = 0
            return np.maximum(abs(lb), abs(ub)) / weight

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
            now, target = np.broadcast_arrays(now, target)
            for vi, vn, vt in zip(np.broadcast_to(ini, now.shape).flat,
                                  now.flat,
                                  np.broadcast_to(target, now.shape).flat):
                strs.append('{:>12s}{: 16e}{: 16e}{: 16e}{: 16e}'.format(
                    name, vi, vn, vt, vn - vt))
        return '\n'.join(strs)


class ElementConstraints(Constraints):
    """Base class for position-related constraints: handle the refpoints
    of each target"""
    def __init__(self, ring, *args, **kwargs):
        self.nelems = len(ring)
        self.refs = []
        self.refpts = bool_refpts([], self.nelems)
        super(ElementConstraints, self).__init__(*args, **kwargs)

    @staticmethod
    def _arrayaccess(index):
        """Access to array elements"""
        if index is None:
            def getv(x):
                return x
        else:
            def getv(x):
                return x[index]
        return getv

    @staticmethod
    def _recordaccess(index):
        """Access to optics parameters"""
        if index is None:
            getf = getattr
        else:
            if isinstance(index, tuple):
                idx = (Ellipsis,)+index
            else:
                idx = (Ellipsis, index)

            def getf(lindata, attrname):
                return getattr(lindata, attrname)[idx]
        return getf

    def add(self, fun, target, refpts=None, **kwargs):
        ref = bool_refpts(refpts, self.nelems)
        # Store the new refpoint
        self.refs.append(ref)
        # Update the union of all refpoints
        self.refpts = np.stack((self.refpts, ref), axis=0).any(axis=0)
        super(ElementConstraints, self).add(fun, target, **kwargs)

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


class LinoptConstraints(ElementConstraints):
    """Container for linear optics constraints:
      - a constraint can be set on any result of at.get_optics
      - constraints are added to the container with the LinoptConstraints.add
        method.

      at.get_optics is called once before the evaluation of all constraints

      Example:
          cnstrs = LinoptConstraints(ring, dp=0.01, coupled=False)

          # Add a beta H (beta[0]) constraint at location ref_inj
          cnstrs.add('beta_x_inj', 'beta', 18.0, refpts=ref_inj, index=0)

          # Add a tune constraint
          cnstrs.add('tunes', 0.44, index=0, weight=0.01)

          # Add a chromaticity constraint (both planes)
          cnstrs.add('chroms', [0.0 0.0])

          # define a constraint of phase advances between 2 points
          def mu_diff(lindata, tune, chrom):
              delta_mu = (lindata[1].mu - lindata[0].mu)/(2*np.pi)
              return delta_mu % 1.0

          # Add a H phase advance constraint, giving the desired locations
          cnstrs.add(mu_diff, 0.5, refpts=[sf0 sf1], index=0)
    """
    def __init__(self, ring, **kwargs):
        """Build a LinoptConstraints container

        KEYWORDS
        dp=0.0          momentum deviation.
        twiss_in=None   Initial twiss parameters for transfer line optics.
                        "lindata" stucture, where only the beta and alpha are
                        required and used.
        orbit=None      Initial trajectory for transfer line
                        ((6,) array)
        method=linopt6  Method used for the analysis of the transfer matrix.  
                        Can be None, at.linopt2, at.linopt4, at.linopt6
                        linopt2:    no longitudinal motion, no H/V coupling,
                        linopt4:    no longitudinal motion, Sagan/Rubin
                                    4D-analysis of coupled motion,
                        linopt6:    with or without longitudinal motion, normal
                                    mode analysis
        """
        self.get_chrom = False
        super(LinoptConstraints, self).__init__(ring, **kwargs)

    def add(self, param, target, refpts=None, index=None, name=None, **kwargs):
        """Add a target to the LinoptConstraints container

        PARAMETERS
            param         2 possibilities:
                          - parameter name: see at.linopt for the name of
                            available parameters. In addition to local optical
                            parameters, 'tunes' and 'chroms' are allowed.
                          - user-supplied parameter evaluation function:
                                value = param(lindata, tune, chrom)
                            lindata contains the optics parameters at all the
                              specified refpoints
                            value is the constrained parameter value
                              (scalar or array).
            target        desired value.

        KEYWORDS
            refpts=None   location of the constraint. Several locations may be
                          given to apply the same constraint at several points.
            index=None    index in the parameter array. If None, the full array
                          is used.
            name=None     name of the constraint. If None, name is generated
                          from param and index.
            weight=1.0    weight factor: the residual is (value-target)/weight.
            bounds=(0,0)  lower and upper bounds. The parameter is constrained
                          in the interval [target-low_bound target+up_bound]
            UseInteger    Match integer part of mu, much slower as the optics
                          calculation is done for all refpts

        The target, weight and bounds values must be broadcastable to the shape
        of value.
        """
        getf = self._recordaccess(index)
        getv = self._arrayaccess(index)
        use_integer = kwargs.pop('UseInteger', False)

        if name is None:                # Generate the constraint name
            name = param.__name__ if callable(param) else param
            if index is not None:
                name = '{}_{}'.format(name, index)

        if callable(param):
            def fun(refdata, tune, chrom):
                return getv(param(refdata, tune, chrom))
            # self.refpts[:] = True     # necessary not to miss 2*pi jumps
            self.get_chrom = True       # fun may use dispersion or chroma
        elif param == 'tunes':
            # noinspection PyUnusedLocal
            def fun(refdata, tune, chrom):
                return getv(tune)
            refpts = []
        elif param == 'chroms':
            # noinspection PyUnusedLocal
            def fun(refdata, tune, chrom):
                return getv(chrom)
            refpts = []
            self.get_chrom = True       # slower but necessary
        else:
            # noinspection PyUnusedLocal
            def fun(refdata, tune, chrom):
                if param == 'mu':
                    return getf(refdata, param)%(2*np.pi)
                elif param == 'mu' and use_integer:
                    self.refpts[:] = True # necessary not to miss 2*pi jumps
                return getf(refdata, param)

        super(LinoptConstraints, self).add(fun, target, refpts, name=name,
                                           **kwargs)

    def compute(self, ring, *args, **kwargs):
        """Optics computation before evaluation of all constraints"""
        ld0, bd, ld = get_optics(ring, refpts=self.refpts,
                                 get_chrom=self.get_chrom, **kwargs)
        return (ld[ref[self.refpts]] for ref in self.refs), \
               (bd.tune[:2], bd.chromaticity)


class OrbitConstraints(ElementConstraints):
    """Container for orbit constraints:
    The closed orbit can be handled with LinoptConstraints, but for problems
    which do not involve parameters other than orbit, like steering or
    orbit bumps, OrbitConstraints is much faster.

      at.find_orbit is called once before the evaluation of all constraints

    Example:
        cnstrs = OrbitConstraints(ring, dp=0.01)

        # Add a bump (x=-0.004, x'=0) constraint at location ref_inj
        cnstrs.add([-0.004, 0.0], refpts=ref_inj, index=slice(2))
    """
    def __init__(self, ring, *args, **kwargs):
        """Build a OrbitConstraints container

        KEYWORDS
        dp=0            Momentum deviation, when radiation is OFF
        dct=0            Path lengthening, when radiation ids OFF
        orbit=None          Initial trajectory for transfer line: (6,) array
        """
        if ring.radiation:
            args.pop('dp',0.0)
            args.pop('dct',0.0)
        super(OrbitConstraints, self).__init__(ring, *args, **kwargs)

    def add(self, target, refpts=None, index=None, name=None, **kwargs):
        """Add a target to the OrbitConstraints container

        PARAMETERS
            target        desired value.

        KEYWORDS
            refpts=None   location of the constraint. Several locations may be
                          given to apply the same constraint at several points.
            index=None    index in the orbit vector. If None, the full orbit
                          is used. Example:
                            index=0         # x
                            index=2         # z
                            index=slice(4)  # x, x', z, z'
            name='orbit'  name of the constraint.
            weight=1.0    weight factor: the residual is (value-target)/weight.
            bounds=(0,0)  lower and upper bounds. The parameter is constrained
                          in the interval [target-low_bound target+up_bound]

        The target, weight and bounds values must be broadcastable to the shape
        of value.
        """

        if name is None:                # Generate the constraint name
            name = 'orbit' if index is None else 'orbit_{}'.format(index)

        fun = self._arrayaccess(index)
        super(OrbitConstraints, self).add(fun, target, refpts, name=name,
                                           **kwargs)

    def compute(self, ring, *args, **kwargs):
        """Orbit computation before evaluation of all constraints"""
        orbit0, orbit = find_orbit(ring, refpts=self.refpts, **kwargs)
        return (orbit[ref[self.refpts]].T for ref in self.refs), ()


class EnvelopeConstraints(ElementConstraints):
    """Container for envelope constraints:
      - a constraint can be set on any result of at.ohmi_envelope
      - constraints are added to the container with the EnvelopeConstraints.add
        method.

      at.ohmi_envelope is called once before the evaluation of all constraints

      Example:
          cnstrs = EnvelopeConstraints(ring)
    """
    def __init__(self, ring):
        """Build a EnvelopeConstraints container"""
        super(EnvelopeConstraints, self).__init__(ring, rad=True)

    def add(self, param, target, refpts=None, index=None, name=None, **kwargs):
        """Add a target to the EnvelopeConstraints container

        PARAMETERS
            param         2 possibilities:
                          - parameter name: see at.ohmi_envelope for the name of
                            available parameters. In addition to local
                            parameters, 'tunes', 'damping_rates',
                            'mode_matrices' and 'mode_emittance' are allowed.
                          - user-supplied parameter evaluation function:
                                value = prm(emit_data, beam_data)
                            emit_data contains the emittance data at all the
                              specified refpoints
                            value is the constrained parameter value
                              (scalar or array).
            target        desired value.

        KEYWORDS
            refpts=None   location of the constraint. Several locations may be
                          given to apply the same constraint at several points.
            index=None    index in the parameter array. If None, the full array
                          is used.
            name=None     name of the constraint. If None, name is generated
                          from param and index.
            weight=1.0    weight factor: the residual is (value-target)/weight.
            bounds=(0,0)  lower and upper bounds. The parameter is constrained
                          in the interval [target-low_bound target+up_bound]

        The target, weight and bounds values must be broadcastable to the shape
        of value.
        """
        def beam_constraint(prm):
            # noinspection PyUnusedLocal
            def beamfunc(refdata, beam_data):
                return getv(beam_data[prm])
            return beamfunc

        getf = self._recordaccess(index)
        getv = self._arrayaccess(index)

        if name is None:                # Generate the constraint name
            name = param.__name__ if callable(param) else param
            if index is not None:
                name = '{}_{}'.format(name, index)

        if callable(param):
            def fun(refdata, beam_data):
                return getv(param(refdata, beam_data))
        elif param in ['tunes', 'damping_rates', 'mode_matrices',
                       'mode_emittances']:
            fun = beam_constraint(param)
            refpts = []
        else:
            # noinspection PyUnusedLocal
            def fun(refdata, beam_data):
                return getf(refdata, param)

        super(EnvelopeConstraints, self).add(fun, target, refpts, name=name,
                                             **kwargs)

    def compute(self, ring, *args, **kwargs):
        """Optics computation before evaluation of all constraints"""
        em0, beamdata, em = ohmi_envelope(ring, refpts=self.refpts, **kwargs)
        return (em[ref[self.refpts]] for ref in self.refs), (beamdata,)


def match(ring, variables, constraints, verbose=2, max_nfev=1000,
          diff_step=1.0e-10, method=None, copy=True):
    """Perform matching of constraints by varying variables

    PARAMETERS
        ring                ring lattice or transfer line
        variables           sequence of Variable objects
        constraints         sequence of Constraints objects

    KEYWORDS
        verbose=2           See scipy.optimize.least_squares
        max_nfev=1000           "
        diff_step=1.0e-10       "
    """
    def fun(vals):
        for value, variable in zip(vals, variables):
            variable.set(ring1, value)

        c = [cons.evaluate(ring1) for cons in constraints]
        return np.concatenate(c, axis=None)

    if copy:
        ring1 = ring.copy()         # Make a shallow copy of ring
        for var in variables:       # Make a deep copy of varying elements
            for ref in uint32_refpts(var.refpts, len(ring1)):
                ring1[ref] = ring1[ref].deepcopy()
    else:
        ring1 = ring  

    aaa = [(var.get(ring1), var.bounds) for var in variables]
    vini, bounds = zip(*aaa)
    bounds = np.array(bounds).T

    cini = [cst.values(ring1) for cst in constraints]
    ntargets = sum(np.size(a) for a in chain.from_iterable(cini))

    if method is None:
        if np.all(abs(bounds) == np.inf) and ntargets >= len(variables):
            method = 'lm'
        else:
            method = 'trf'
    if verbose >= 1:
        print('\n{} constraints, {} variables, using method {}\n'.
              format(ntargets, len(variables), method))

    least_squares(fun, vini, bounds=bounds, verbose=verbose, max_nfev=max_nfev,
                  method=method, diff_step=diff_step)

    if verbose >= 1:
        print(Constraints.header())
        for cst, ini in zip(constraints, cini):
            print(cst.status(ring1, initial=ini))


        print(Variable.header())
        for var, vini in zip(variables, vini):
            print(var.status(ring1, vini=vini))

    return ring1
