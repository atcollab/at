"""
Classes for matching variables and constraints (obsolete)

**These classes and functions are now obsolete** and should be replaced by
:doc:`variables <at.lattice.lattice_variables>`,
:doc:`observables <at.latticetools.observables>`,
:doc:`matching <at.latticetools.matching>`.
"""
from __future__ import annotations
from itertools import chain
import numpy as np
from collections.abc import Sequence, Callable
from typing import Optional, Union
from scipy.optimize import least_squares
from itertools import repeat
from at.lattice import Lattice, Refpts, bool_refpts
from at.physics import get_optics, ohmi_envelope, find_orbit

__all__ = ['Variable', 'ElementVariable', 'Constraints', 'ElementConstraints',
           'LinoptConstraints', 'OrbitConstraints', 'EnvelopeConstraints',
           'match']


class Variable(object):
    """A :py:class:`Variable` is a scalar value acting on a lattice through the
    user-defined functions setfun and getfun

    Parameters:
        setfun:     User-defined function for setting the Variable. Called as:

          :code:`setfun(ring, value, *args, **kwargs)`

          where :code:`value` is the scalar value to apply

          The positional and keyword parameters come from the
          :py:class:`Variable` initialisation
        getfun:     User-defined function for retrieving the actual value of
          the variable: Called as:

          :code:`value = getfun(ring, *args, **kwargs)`

          The positional and keyword parameters come from the
          :py:class:`Variable` initialisation
        name:       Name of the Variable; Default: ``''``
        bounds:     Lower and upper bounds of the variable value
        fun_args:       Positional arguments transmitted to ``setfun`` and
          ``getfun`` functions

    Keyword Args:
        **fun_kwargs:   Keyword arguments transmitted to ``setfun``and
          ``getfun`` functions
    """
    def __init__(self, setfun: Callable, getfun: Callable,
                 name: str = '',
                 bounds: tuple[float, float] = (-np.inf, np.inf),
                 fun_args: tuple = (), **fun_kwargs):
        self.setfun = setfun
        self.getfun = getfun
        self.name = name
        self.bounds = bounds
        self.args = fun_args
        self.kwargs = fun_kwargs
        super(Variable, self).__init__()

    def set(self, ring: Lattice, value):
        """Set the variable value"""
        self.setfun(ring, value, *self.args, **self.kwargs)

    def get(self, ring: Lattice):
        """Get the actual variable value"""
        return self.getfun(ring, *self.args, **self.kwargs)

    @staticmethod
    def header():
        return '\n{:>12s}{:>13s}{:>16s}{:>16s}\n'.format(
            'Name', 'Initial', 'Final ', 'Variation')

    def status(self, ring: Lattice, vini=np.nan):
        vnow = self.get(ring)
        return '{:>12s}{: 16e}{: 16e}{: 16e}'.format(
            self.name, vini, vnow, (vnow - vini))


class ElementVariable(Variable):
    """An :py:class:`ElementVariable` is:

    * a scalar attribute or
    * an element of an array attribute

    of one or several :py:class:`.Element` (s) of a lattice.

    Parameters:
        refpts:     Location of variable :py:class:`.Element` (s)
        attname:    Attribute name
        index:      Index in the attribute array. Use :py:obj:`None` for
          scalar attributes
        name:       Name of the Variable; Default: ``''``
        bounds:     Lower and upper bounds of the variable value
    """

    def __init__(self, refpts: Refpts, attname: str,
                 index: Optional[int] = None,
                 name: str = '',
                 bounds: tuple[float, float] = (-np.inf, np.inf)):
        setf, getf = self._access(index)

        def setfun(ring, value):
            for elem in ring.select(refpts):
                setf(elem, attname, value)

        def getfun(ring):
            values = np.array([getf(elem, attname) for elem in
                               ring.select(refpts)])
            return np.average(values)

        super(ElementVariable, self).__init__(setfun, getfun, name=name,
                                              bounds=bounds)
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

    * a constraint is defined by a user-defined evaluation function,
    * constraints are added to the container with the :py:meth:`add` method.

    Parameters:
        *args:      Positional arguments sent to the evaluation functions
          of all the embedded constraints

    Keyword Args:
        **kwargs:   Keyword arguments sent to the evaluation functions
          of all the embedded constraints

    Examples:
        Define an evaluation function for the ring circumference:

        >>> def circ_fun(ring):
        ...     return ring.get_s_pos(len(ring) + 1)

        Define an evaluation function for the momentum compaction factor:

        >>> def mcf_fun(ring):
        ...     return ring.get_mcf()

        Construct the container:

        >>> cnstrs = Constraints()

        Add the two constraints:

        >>> cnstrs.add(circ_fun, 850.0)
        >>> cnstrs.add(mcf_fun, 1.0e-4, weight=0.1)
    """
    def __init__(self, *args, **kwargs):
        self.name = []
        self.fun = []
        self.target = []
        self.weight = []
        self.lbound = []
        self.ubound = []
        self.rad = kwargs.pop('rad', False)
        self.args = args
        self.kwargs = kwargs

    def add(self, fun: Callable, target, name: Optional[str] = None,
            weight=1.0, bounds=(0.0, 0.0)):
        """Add a target to the :py:class:`Constraints` container

        .. highlight:: python

        Parameters:
            fun:          Evaluation function. Called as:

              :code:`value = fun(ring, *args, **kwargs)`

              ``value`` is the constrained parameter value,

              ``value`` may be a scalar or an array.

              The positional and keyword parameters come from
              the :py:class:`Constraints` initialisation
            target:       Desired value.
            name:         Name of the constraint. If :py:obj:`None`, a ``name``
              is generated from the name of the evaluation function
            weight:       Weight factor: the residual is
              :code:`(value-target)/weight`
            bounds:       Lower and upper bounds. The parameter is constrained
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

    def values(self, ring: Lattice):
        """Return the list of actual parameter values"""
        return [fun(ring, *self.args, **self.kwargs) for fun in self.fun]

    def evaluate(self, ring: Lattice):
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

    def status(self, ring: Lattice, initial=None):
        """Return a string giving the actual state of constraints"""
        if initial is None:
            initial = repeat(np.nan)
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
    def __init__(self, ring: Lattice, *args, **kwargs):
        self.nelems = len(ring)
        self.refs = []
        self.refpts = ring.bool_refpts([])
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

    def add(self, fun: Callable, target,
            refpts: Optional[Refpts] = None, **kwargs):
        """Add a target to the :py:class:`ElementConstraints` container

        Parameters:
            fun:          Evaluation function. Called as:

              :code:`value = fun(ring, *args, **kwargs)`

              ``value`` is the constrained parameter value

              ``value`` may be a scalar or an array.

              the positional and keyword parameters come from
              the :py:class:`ElementConstraints` initialisation
            target:       Desired value.
            refpts:       Location of the constraint

        Keyword Args:
            name:         Name of the constraint. If :py:obj:`None`, a ``name``
              is generated from the name of the evaluation function
            weight:       Weight factor: the residual is
              :code:`(value-target)/weight`
            bounds:       Lower and upper bounds. The parameter is constrained
              in the interval [target-low_bound target+up_bound]

        The "target", "weight" and "bounds" input must be broadcastable to the
        shape of "value".
        """
        ref = bool_refpts(refpts, self.nelems)
        # Store the new refpoint
        self.refs.append(ref)
        # Update the union of all refpoints
        self.refpts = np.stack((self.refpts, ref), axis=0).any(axis=0)
        super(ElementConstraints, self).add(fun, target, **kwargs)

    def values(self, ring: Lattice):
        # Single optics computation
        vloc, vglob = self.compute(ring, *self.args, **self.kwargs)
        # Evaluate all constraints
        return [fun(loc, *vglob) for fun, loc in zip(self.fun, vloc)]

    def compute(self, ring: Lattice, *args, **kwargs):
        """Dummy computation. Compute must return:
        - an iterator over local data for each target
        - a tuple of global data
        """
        return repeat(None), ()


class LinoptConstraints(ElementConstraints):
    # noinspection PyUnresolvedReferences
    """Container for linear optics constraints:

    * a constraint can be set on any result of at.get_optics
    * constraints are added to the container with the LinoptConstraints.add
      method.

    :py:func:`.get_optics` is called once before the evaluation of all
    constraints

    Parameters:
        ring:       Lattice description

    Keyword Args:
        dp (float):   Momentum deviation.
        dct (float):  Path lengthening. If specified, ``dp`` is
          ignored and the off-momentum is deduced from the path lengthening.
        orbit (Optional[Orbit]): Avoids looking for the closed orbit if is
          already known ((6,) array)
        twiss_in:   Initial twiss parameters for transfer line optics.
          See :py:func:`.linopt6`
        method (Callable):  Method used for the analysis of the
          transfer matrix. Can be :py:obj:`~.linear.linopt2`,
          :py:obj:`~.linear.linopt4`, :py:obj:`~.linear.linopt6`

          * :py:obj:`~.linear.linopt2`: No longitudinal motion,
                                        no H/V coupling,
          * :py:obj:`~.linear.linopt4`: No longitudinal motion, Sagan/Rubin
            4D-analysis of coupled motion,
          * :py:obj:`~.linear.linopt6` (default): With or without longitudinal
            motion, normal mode analysis

    Example:

        >>> cnstrs = LinoptConstraints(ring, dp=0.01, coupled=False)

        Add a beta x (beta[0]) constraint at location ref_inj:

        >>> cnstrs.add('beta', 18.0, refpts=ref_inj,
                       name='beta_x_inj', index=0)

        Add an horizontal tune (tunes[0]) constraint:

        >>> cnstrs.add('tunes', 0.44, index=0, weight=0.01)

        Add a chromaticity constraint (both planes):

        >>> cnstrs.add('chroms', [0.0 0.0])

        Define a constraint of phase advances between 2 points:

        >>> def mu_diff(lindata, tune, chrom):
        ...     delta_mu = (lindata[1].mu - lindata[0].mu)/(2*np.pi)
        ...     return delta_mu % 1.0

        Add a H phase advance constraint, giving the desired locations:

        >>> cnstrs.add(mu_diff, 0.5, refpts=[sf0 sf1], index=0)
        """
    def __init__(self, ring: Lattice, **kwargs):
        self.get_chrom = False
        super(LinoptConstraints, self).__init__(ring, **kwargs)

    def add(self, param, target, refpts: Optional[Refpts] = None,
            index: Optional[Union[int, slice]] = None,
            name: Optional[str] = None, **kwargs):
        """Add a target to the LinoptConstraints container

        Parameters:
            param:         2 possibilities:

              * parameter name: see :py:func:`.linopt6` for the
                name of available parameters. In addition to local
                optical parameters, ``'tunes'`` and ``'chroms'``
                are allowed.
              * user-supplied parameter evaluation function:

                :code:`value = param(lindata, tune, chrom)`

                ``lindata`` contains the optics parameters at all the
                specified refpoints
                ``value`` is the constrained parameter value
                (scalar or array).
            target:       desired value.
            refpts:       location of the constraint. Several locations may be
              given to apply the same constraint at several points.
            index:        index in the parameter array. If :py:obj:`None`,
              the full array is used.
            name:         name of the constraint. If :py:obj:`None`, name is
              generated from ``param`` and ``index``.

        Keyword Args:
            weight:       Weight factor: the residual is
              :code:`(value-target)/weight`
            bounds:       lower and upper bounds. The parameter is constrained
                          in the interval [target-low_bound target+up_bound]
            UseInteger:   Match integer part of mu, much slower as the optics
                          calculation is done for all refpts

        The target, weight and bounds values must be broadcastable to the shape
        of value.
        """
        getf = self._recordaccess(index)
        getv = self._arrayaccess(index)
        use_integer = kwargs.pop('UseInteger', False)
        norm_mu = {'mu': 1, 'mun': 2*np.pi}

        if name is None:                # Generate the constraint name
            name = param.__name__ if callable(param) else param
            if index is not None:
                name = '{}_{}'.format(name, index)

        if callable(param):
            def fun(refdata, tune, chrom):
                return getv(param(refdata, tune, chrom))
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
        elif param == 'mu' or param == 'mun':
            # noinspection PyUnusedLocal
            def fun(refdata, tune, chrom):
                if use_integer:
                    return getf(refdata, 'mu') / norm_mu[param]
                else:
                    return (getf(refdata, 'mu') % (2*np.pi)) / norm_mu[param]
            if use_integer:
                self.refpts[:] = True  # necessary not to miss 2*pi jumps
            else:
                target = target % (2 * np.pi / norm_mu[param])
        else:
            # noinspection PyUnusedLocal
            def fun(refdata, tune, chrom):
                return getf(refdata, param)

        super(LinoptConstraints, self).add(fun, target, refpts, name=name,
                                           **kwargs)

    def compute(self, ring: Lattice, *args, **kwargs):
        """Optics computation before evaluation of all constraints"""
        ld0, bd, ld = get_optics(ring, refpts=self.refpts,
                                 get_chrom=self.get_chrom, **kwargs)
        return (ld[ref[self.refpts]] for ref in self.refs), \
               (bd.tune[:2], bd.chromaticity)


class OrbitConstraints(ElementConstraints):
    # noinspection PyUnresolvedReferences
    """Container for orbit constraints:
    The closed orbit can be handled with :py:class`LinoptConstraints`, but for
    problems which do not involve parameters other than orbit, like steering or
    orbit bumps, :py:class:`OrbitConstraints` is much faster.

    :py:func:`.find_orbit` is called once before the evaluation of all
    constraints

    Parameters:
        ring:       Lattice description

    Keyword Args:
        dp (float):   Momentum deviation.
        dct (float):  Path lengthening. If specified, ``dp`` is
          ignored and the off-momentum is deduced from the path lengthening.
        orbit (Optional[Orbit]): Avoids looking for the closed orbit if is
          already known ((6,) array)

    Example:
        >>> cnstrs = OrbitConstraints(ring, dp=0.01)

        Add a bump (x=-0.004, x'=0) constraint at location ref_inj

        >>> cnstrs.add([-0.004, 0.0], refpts=ref_inj, index=slice(2))
    """
    def __init__(self, ring: Lattice, *args, **kwargs):
        if ring.radiation:
            kwargs.pop('dp', 0.0)
            kwargs.pop('dct', 0.0)
        super(OrbitConstraints, self).__init__(ring, *args, **kwargs)

    def add(self, target, refpts: Optional[Refpts] = None,
            index: Optional[Union[int, slice]] = None,
            name: Optional[str] = None, **kwargs):
        """Add a target to the OrbitConstraints container

        Parameters:
            target:       desired value.
            refpts:       location of the constraint. Several locations may be
              given to apply the same constraint at several points.
            index:        index in the orbit vector. If :py:obj:`None`, the
              full orbit is used.
            name:         name of the constraint. Default: ``'orbit'``

        Keyword Args:
            weight:       Weight factor: the residual is
              :code:`(value-target)/weight`
            bounds:       lower and upper bounds. The parameter is constrained
                          in the interval [target-low_bound target+up_bound]

        The target, weight and bounds values must be broadcastable to the shape
        of value.
        """

        if name is None:                # Generate the constraint name
            name = 'orbit' if index is None else 'orbit_{}'.format(index)

        fun = self._arrayaccess(index)
        super(OrbitConstraints, self).add(fun, target, refpts, name=name,
                                          **kwargs)

    def compute(self, ring: Lattice, *args, **kwargs):
        """Orbit computation before evaluation of all constraints"""
        orbit0, orbit = find_orbit(ring, refpts=self.refpts, **kwargs)
        return (orbit[ref[self.refpts]].T for ref in self.refs), ()


class EnvelopeConstraints(ElementConstraints):
    """Container for envelope constraints:

    * a constraint can be set on any result of :py:func:`.ohmi_envelope`,
    * constraints are added to the container with the :py:meth:`add` method.

    :py:func:`.ohmi_envelope` is called once before the evaluation of all
    constraints.

    Parameters:
        ring:       Lattice description
    """
    def __init__(self, ring: Lattice):
        super(EnvelopeConstraints, self).__init__(ring, rad=True)

    def add(self, param, target, refpts: Optional[Refpts] = None,
            index: Optional[Union[int, slice]] = None,
            name: Optional[str] = None, **kwargs):
        """Add a target to the :py:class:`EnvelopeConstraints` container

        Parameters:
            param:        2 possibilities:

              - parameter name: see :py:func:`.ohmi_envelope` for the
                name of available parameters. In addition to local parameters,
                ``'tunes'``, ``'damping_rates'``, ``'mode_matrices'`` and
                ``'mode_emittance'`` are allowed.
              - user-supplied parameter evaluation function:

                :code:`value = param(emit_data, beam_data)`

                ``emit_data`` contains the emittance data at all the
                specified refpoints
                ``value`` is the constrained parameter value
                (scalar or array).
            target:       desired value.
            refpts:       location of the constraint. Several locations may be
              given to apply the same constraint at several points.
            index:        index in the parameter array. If :py:obj:`None`,
              the full array is used.
            name:         name of the constraint. If :py:obj:`None`, name is
              generated from ``param`` and ``index``.

        Keyword Args:
            weight:       Weight factor: the residual is
              :code:`(value-target)/weight`
            bounds:       lower and upper bounds. The parameter is constrained
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

    def compute(self, ring: Lattice, *args, **kwargs):
        """Optics computation before evaluation of all constraints"""
        em0, beamdata, em = ohmi_envelope(ring, refpts=self.refpts, **kwargs)
        return (em[ref[self.refpts]] for ref in self.refs), (beamdata,)


def match(ring: Lattice, variables: Sequence[Variable],
          constraints: Sequence[Constraints], verbose: int = 2,
          max_nfev: int = 1000,
          diff_step: float = 1.0e-10,
          method=None, copy: bool = True):
    """Perform matching of constraints by varying variables (obsolete)

    Parameters:
        ring:               Lattice description
        variables:          sequence of Variable objects
        constraints:        sequence of Constraints objects
        verbose:            Print additional information
        max_nfev:           Maximum number of evaluations
        diff_step:          Convergence threshold
        method:
        copy:
    """
    def fun(vals):
        for value, variable in zip(vals, variables):
            variable.set(ring1, value)

        c = [cons.evaluate(ring1) for cons in constraints]
        return np.concatenate(c, axis=None)

    if copy:
        # Make a shallow copy of ring
        ring1 = ring.copy()
        varpts = ring.bool_refpts([])
        # build the list of variable elements
        for var in variables:
            if isinstance(var, ElementVariable):
                varpts |= ring.bool_refpts(var.refpts)
        # make a deep copy of all the variable elements
        for ref in ring.uint32_refpts(varpts):
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
