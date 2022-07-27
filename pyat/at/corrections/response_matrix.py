"""
Classes to compute arbitrary response matrices
"""

import numpy
from at.matching import Variable, ElementVariable


class Observables(object):
    """Container for generic observable:

    * a observable is defined by a user-defined evaluation function,
    * observable are added to the container with the :py:meth:`add` method.

    Parameters:
        *args:      Positional arguments sent to the evaluation functions
          of all the embedded constraints

    Keyword Args:
        **kwargs:   Keyword arguments sent to the evaluation functions
          of all the embedded constraints

    Examples:
        Define an evaluation function for the horizontal orbit:

        >>> def orbith(ring, refpts=None):
        ...     return ring.find_orbit(refpts=refpts)[0,:]

        Define an evaluation function for the vertical orbit:

        >>> def orbitv(ring, refpts=None):
        ...     return ring.find_orbit(refpts=refpts)[2,:]

        Construct the container:

        >>> obs = Observables(refpts=refpts)

        Add the two constraints:

        >>> obs.add(orbith, "Hor. orbit")
        >>> obs.add(orbitv, "Ver. orbit")
    """
    def __init__(self, *args, **kwargs):
        self.name = []
        self.fun = []
        self.args = args
        self.kwargs = kwargs

    def add(self, fun: Callable, name: Optional[str] = None):
        """Add an observable to the :py:class:`Observables` container

        .. highlight:: python

        Parameters:
            fun:          Evaluation function. Called as:

              :code:`value = fun(ring, *args, **kwargs)`

              ``value`` is the constrained parameter value,

              ``value`` may be a scalar or an array.

              The positional and keyword parameters come from
              the :py:class:`Observables` initialisation
            name:         Name of the observable. If :py:obj:`None`, a ``name``
              is generated from the name of the evaluation function
        """
        if name is None:                # Generate the constraint name
            name = fun.__name__
        self.name.append(name)
        self.fun.append(fun)

    def values(self, ring: Lattice):
        """Return the list of actual parameter values"""
        return [fun(ring, *self.args, **self.kwargs) for fun in self.fun]
        

class ElementObservables(Observables):
    """Base class for position-related constraints: handle the refpoints
    of each target"""
    def __init__(self, ring: Lattice, *args, **kwargs):
        self.nelems = len(ring)
        self.refs = []
        self.refpts = ring.bool_refpts([])
        super(ElementObservables, self).__init__(*args, **kwargs)

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
        """Add an observable to the :py:class:`ElementObservables` container

        Parameters:
            fun:          Evaluation function. Called as:

              :code:`value = fun(ring, *args, **kwargs)`

              ``value`` is the constrained parameter value

              ``value`` may be a scalar or an array.

              the positional and keyword parameters come from
              the :py:class:`ElementObservables` initialisation
            refpts:       Location of the observable

        Keyword Args:
            name:         Name of the observable. If :py:obj:`None`, a ``name``
              is generated from the name of the evaluation function
        """
        ref = bool_refpts(refpts, self.nelems)
        # Store the new refpoint
        self.refs.append(ref)
        # Update the union of all refpoints
        self.refpts = np.stack((self.refpts, ref), axis=0).any(axis=0)
        super(ElementObservable, self).add(fun, target, **kwargs)

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
        
class OrbitObservables(ElementObservables):
    # noinspection PyUnresolvedReferences
    """Container for orbit observables:
    The closed orbit can be handled with :py:class`LinoptObservables`, but for
    problems which do not involve parameters other than orbit, like steering or
    orbit bumps, :py:class:`OrbitObservables` is much faster.

    Parameters:
        ring:       Lattice description

    Keyword Args:
        dp (float):   Momentum deviation.
        dct (float):  Path lengthening. If specified, ``dp`` is
          ignored and the off-momentum is deduced from the path lengthening.

    Example:
        >>> obss = OrbitObservables(ring, dp=0.01)

        Add an observable at location ref_inj

        >>> obs.add(refpts=ref_inj, index=slice(2))
    """
    def __init__(self, ring: Lattice, *args, **kwargs):
        if ring.radiation:
            kwargs.pop('dp', 0.0)
            kwargs.pop('dct', 0.0)
        super(OrbitObservables, self).__init__(ring, *args, **kwargs)

    def add(self, refpts: Optional[Refpts] = None,
            index: Optional[Union[int, slice]] = None,
            name: Optional[str] = None):
        """Add an observable to the OrbitObservables container

        Parameters:
            refpts:       location of the constraint. Several locations may be
              given to apply the same observable at several points.
            index:        index in the orbit vector. If :py:obj:`None`, the
              full orbit is used.
            name:         name of the observable. Default: ``'orbit'``
        """

        if name is None:                # Generate the constraint name
            name = 'orbit' if index is None else 'orbit_{}'.format(index)

        fun = self._arrayaccess(index)
        super(OrbitObservables, self).add(fun, refpts, name=name)

    def compute(self, ring: Lattice, *args, **kwargs):
        """Orbit computation before evaluation of all observables"""
        orbit0, orbit = find_orbit(ring, refpts=self.refpts, **kwargs)
        return (orbit[ref[self.refpts]].T for ref in self.refs), ()
