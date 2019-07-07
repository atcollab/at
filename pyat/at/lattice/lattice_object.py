"""Lattice object

The methods implemented in this module are internal to the 'lattice' package.
This is necessary to ensure that the 'lattice' package is independent
from other AT packages.

Other Lattice methods are implemented in other AT packages and are available
as soon as the package is imported. The 'tracking' and 'physics' packages are
automatically imported.

As an example, see the at.physics.orbit module
"""
import sys
import copy
import numpy
import math
import itertools
from warnings import warn
from at.lattice import AtError, AtWarning
from at.lattice import elements, get_s_pos, get_elements, uint32_refpts

__all__ = ['Lattice']


class Lattice(list):
    """Lattice object
    An AT lattice is a sequence of AT elements.
    A Lattice accepts extended indexing (as a numpy ndarray).

    ATTRIBUTES
        name        Name of the lattice
        energy      Particle energy
        periodicity Number of super-periods to describe the full ring

    To reduce the inter-package dependencies, some methods of the
    lattice object are defined in other AT packages, in the module where
    the underlying function is implemented.
    """
    _1st_attributes = ('name', 'energy', 'periodicity')

    def __init__(self, elems=None, name=None, energy=None, periodicity=None,
                 **kwargs):
        """Create a new lattice object

        INPUT
            elems:          any iterable of AT elements

        KEYWORDS
            name            Name of the lattice
            energy          Energy of the lattice
            periodicity     Number of periods

            all other keywords will be set as Lattice attributes
        """
        # noinspection PyShadowingNames
        def check_elems(elems, attrs):
            """Check the radiation status of a lattice"""
            radiate = False
            for idx, elem in enumerate(elems):
                if isinstance(elem, elements.Element):
                    if (elem.PassMethod.endswith('RadPass') or
                            elem.PassMethod.endswith('CavityPass')):
                        radiate = True
                    yield elem
                else:
                    warn(AtWarning('item {0} ({1}) is not an AT element: '
                                   'ignored'.format(idx, elem)))
            attrs['_radiation'] = radiate

        if elems is None:
            elems = []

        if name is not None:
            kwargs['name'] = name
        if energy is not None:
            kwargs['energy'] = energy
        if periodicity is not None:
            kwargs['periodicity'] = periodicity

        if isinstance(elems, Lattice):
            attrs = vars(elems).copy()
            attrs.update(kwargs)
        else:
            attrs = kwargs

        # reset the range of interest
        attrs.pop('_s_range', None)
        attrs.pop('_i_range', None)
        # set default values
        attrs.setdefault('name', '')
        attrs.setdefault('periodicity', 1)
        if 'energy' not in attrs:
            raise AtError('Lattice energy is not defined')
        if '_radiation' not in attrs:
            elems = check_elems(elems, attrs)

        super(Lattice, self).__init__(elems)

        for key, value in attrs.items():
            setattr(self, key, value)

    def __getitem__(self, key):
        try:
            return super(Lattice, self).__getitem__(key.__index__())
        except (AttributeError, TypeError):
            return Lattice(self.iterator(key), **vars(self))

    if sys.version_info < (3, 0):
        # This won't be defined if version is at least 3.0
        # Called for slices with step != 1
        # noinspection PyTypeChecker
        def __getslice__(self, i, j):
            return self.__getitem__(slice(i, j))

    def __setitem__(self, key, values):
        try:
            super(Lattice, self).__setitem__(key, values)
        except TypeError:
            key = uint32_refpts(key, len(self))
            for i, v in zip(*numpy.broadcast_arrays(key, values)):
                super(Lattice, self).__setitem__(i, v)

    def __delitem__(self, key):
        try:
            super(Lattice, self).__delitem__(key)
        except TypeError:
            key = uint32_refpts(key, len(self))
            for i in reversed(key):
                super(Lattice, self).__delitem__(i)

    def __repr__(self):
        attrs = vars(self).copy()
        k1 = [(k, attrs.pop(k))for k in Lattice._1st_attributes]
        k2 = [(k, v) for k, v in attrs.items() if not k.startswith('_')]
        keys = ', '.join('{0}={1!r}'.format(k, v) for k, v in (k1+k2))
        return 'Lattice({0}, {1})'.format(super(Lattice, self).__repr__(), keys)

    def __str__(self):
        attrs = vars(self).copy()
        k1 = [(k, attrs.pop(k))for k in Lattice._1st_attributes]
        k2 = [(k, v) for k, v in attrs.items() if not k.startswith('_')]
        keys = ', '.join('{0}={1!r}'.format(k, v) for k, v in (k1+k2))
        return 'Lattice(<{0} elements>, {1})'.format(len(self), keys)

    def __add__(self, elems):
        return Lattice(itertools.chain(self, elems), **vars(self))

    def __mul__(self, times):
        return Lattice(itertools.chain(*itertools.repeat(self, times)),
                       **vars(self))

    def iterator(self, key):
        """Iterates over the indices selected by a slice or an array"""
        if isinstance(key, slice):
            indices = key.indices(len(self))
            rg = range(*indices)
        else:
            rg = uint32_refpts(key, len(self))
        return (super(Lattice, self).__getitem__(i) for i in rg)

    def copy(self):
        """Return a shallow copy"""
        return Lattice(self)

    def deepcopy(self):
        """Return a deep copy"""
        return copy.deepcopy(self)

    def slice(self, size=None, slices=1):
        """Create a new lattice by slicing the range of interest into small
        elements

        KEYWORDS
            size=None       Length of a slice. Default: computed from the
                            range and number of points:
                                    sx = (s_max-s_min)/slices.
            slices=1        Number of slices in the specified range. Ignored if
                            size is specified. Default: no slicing

        RETURN
            New Lattice object
       """

        def slice_iter(ibeg, iend):
            for elem in self[:ibeg]:
                yield elem
            for elem in self[ibeg:iend]:
                nslices = int(math.ceil(elem.Length / size))
                if nslices > 1:
                    frac = numpy.ones(nslices) / nslices
                    for el in elem.divide(frac):
                        yield el
                else:
                    yield elem
            for elem in self[iend:]:
                yield elem

        s_range = self.s_range
        if size is None:
            smin, smax = s_range
            size = (smax - smin) / slices

        i1 = self._i_range[0]
        i2 = self._i_range[-1]
        return Lattice(slice_iter(i1, i2), s_range=s_range, **vars(self))

    @property
    def s_range(self):
        """Range of interest. 'None' means the full cell."""
        try:
            return self._s_range
        except AttributeError:
            self.s_range = None
            return self._s_range

    # noinspection PyAttributeOutsideInit
    @s_range.setter
    def s_range(self, value):
        spos = self.get_s_pos(range(len(self) + 1))
        if value is None:
            value = (0.0, spos[-1])
        ok = numpy.flatnonzero(
            numpy.logical_and(spos > value[0], spos < value[1]))
        if len(ok) > 0:
            i1 = max(ok[0] - 1, 0)
            i2 = min(ok[-1] + 1, len(self))
            self._i_range = range(i1, i2 + 1)
        else:
            self._i_range = range(0, 0)
        self._s_range = value

    @property
    def i_range(self):
        """Range of elements inside the range of interest"""
        try:
            i_range = self._i_range
        except AttributeError:
            self.s_range = None
            i_range = self._i_range
        return uint32_refpts(i_range, len(self))

    @property
    def circumference(self):
        """Ring circumference"""
        return self.periodicity * self.get_s_pos(len(self))[0]

    @property
    def voltage(self):
        """Total accelerating voltage"""
        volts = [elem.Voltage for elem in self if
                 isinstance(elem, elements.RFCavity)]
        return self.periodicity * sum(volts)

    @property
    def harmonic_number(self):
        """Harmonic number"""
        harms = [elem.HarmNumber for elem in self if
                 isinstance(elem, elements.RFCavity)]
        return self.periodicity * harms[0] if harms else None

    @property
    def radiation(self):
        """If True, at least one element modifies the beam energy"""
        return self._radiation

    def _radiation_switch(self, cavity_func=None, dipole_func=None,
                          quadrupole_func=None):

        def switch(elfunc, elcheck):
            if elfunc is not None:
                for el in (e for e in self if elcheck(e)):
                    elfunc(el)

        def isdipole(elem):
            return isinstance(elem, elements.Dipole) and (
                    elem.BendingAngle != 0.0)

        switch(cavity_func, lambda e: isinstance(e, elements.RFCavity))
        switch(dipole_func, isdipole)
        switch(quadrupole_func, lambda e: isinstance(e, elements.Quadrupole))

    def radiation_on(self, cavity_pass='CavityPass', dipole_pass='auto',
                     quadrupole_pass=None):
        """
        Turn acceleration and radiation on

        KEYWORDS
            cavity_pass='CavityPass'    PassMethod set on cavities
            dipole_pass='auto'          PassMethod set on dipoles
            quadrupole_pass=None        PassMethod set on quadrupoles

            For PassMethod names, the convention is:
                None            no change
                'auto'          replace *Pass by *RadPass
                anything else   set as it is

        """

        def repfunc(pass_method):
            if pass_method is None:
                ff = None
            elif pass_method == 'auto':
                def ff(elem):
                    if not elem.PassMethod.endswith('RadPass'):
                        elem.PassMethod = ''.join(
                            (elem.PassMethod[:-4], 'RadPass'))
                        elem.Energy = self.energy
            else:
                def ff(elem):
                    elem.PassMethod = pass_method
                    elem.Energy = self.energy
            return ff

        self._radiation_switch(repfunc(cavity_pass), repfunc(dipole_pass),
                               repfunc(quadrupole_pass))
        # noinspection PyAttributeOutsideInit
        self._radiation = True

    def radiation_off(self, cavity_pass='IdentityPass', dipole_pass='auto',
                      quadrupole_pass=None):
        """
        Turn acceleration and radiation off

        KEYWORDS
            cavity_pass='IdentityPass'  PassMethod set on cavities
            dipole_pass='auto'          PassMethod set on dipoles
            quadrupole_pass=None        PassMethod set on quadrupoles

            For PassMethod names, the convention is:
                None            no change
                'auto'          replace *RadPass by *Pass
                anything else   set as it is

        """

        def repfunc(pass_method):
            if pass_method is None:
                ff = None
            elif pass_method == 'auto':
                def ff(elem):
                    if elem.PassMethod.endswith('RadPass'):
                        elem.PassMethod = ''.join(
                            (elem.PassMethod[:-7], 'Pass'))
            else:
                def ff(elem):
                    elem.PassMethod = pass_method
            return ff

        self._radiation_switch(repfunc(cavity_pass), repfunc(dipole_pass),
                               repfunc(quadrupole_pass))
        # noinspection PyAttributeOutsideInit
        self._radiation = False


Lattice.get_elements = get_elements
Lattice.get_s_pos = get_s_pos
