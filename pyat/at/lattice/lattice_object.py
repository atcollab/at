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

TWO_PI_ERROR = 1.E-4

__all__ = ['Lattice', 'type_filter', 'params_filter', 'lattice_filter',
           'no_filter']


class Lattice(list):
    """Lattice object
    An AT lattice is a sequence of AT elements.
    A Lattice accepts extended indexing (as a numpy ndarray).

    Lattice attributes;
        name        Name of the lattice
        energy      Particle energy
        periodicity Number of super-periods to describe the full ring

    Lattice(elems, **params)     Create a new lattice object

        INPUT
            elems:          any iterable of AT elements

        KEYWORDS
            iterator=None   Custom iterator (see below)
            scan=False      Scan elements, looking for energy and periodicity
            name            Name of the lattice
            energy          Energy of the lattice
            periodicity     Number of periods
            *               all other keywords will be set as Lattice attributes

    To reduce the inter-package dependencies, some methods of the
    lattice object are defined in other AT packages, in the module where
    the underlying function is implemented.

    Custom iterators:
    Instead of running through 'elems', the Lattice constructor can use one
    or several custom iterators.

    Lattice(*args, iterator=it, **params)

        The iterator "it" is called as "it(params, *args)" and must return an
            iterator over AT elements for building the lattice. It must also
            fill the "params" dictionary used to set the Lattice attributes.

            params is the dictionary of lattice parameters. It is initialised
                   with the keywords of the lattice constructor. It can be
                   modified by the custom iterator to add or remove parameters.
                   Finally, all remaining parameters will be set as Lattice
                   attributes.
            *args  all positional arguments of the Lattice constructor are sent
                   to the custom iterator.

        An iterator can be:
            - a "generator" which yields elements from scratch.
              Examples: a list, or a file iterator,
            - a "filter" which runs through an input iterator, processes each
              element, possibly adds parameters to the params dictionary
              and yields the processed elements.

        Example of chaining iterators (taken from "load_mat"):

        Lattice(ringparam_filter, matfile_generator, filename
                iterator=params_filter, **params)

        matfile_generator(params, filename)
            opens filename and generates AT elements for each cell of the
            Matlab cell array representing the lattice,

        ringparam_filter(params, matfile_generator, *args)
            runs through matfile_generator(params, *args), looks for RingParam
            elements, fills params with their information and discards them,

        params_filter(params, ringparam_filter, *args)
            runs through ringparam_filter(params, *args), looks for energy and
            periodicity if not yet defined.
    """
    _1st_attributes = ('name', 'energy', 'periodicity')

    def __init__(self, *args, **kwargs):
        """Lattice constructor"""

        iterator = kwargs.pop('iterator', None)
        scan = kwargs.pop('scan', False)
        for key in list(k for k in kwargs.keys() if k.startswith('_')):
            kwargs.pop(key)
        if iterator is None:
            elem_list, = args or [[]]  # accept 0 or 1 argument
            if isinstance(elem_list, Lattice):
                elems = lattice_filter(kwargs, elem_list)
            elif scan:
                elems = params_filter(kwargs, type_filter, elem_list)
            else:
                elems = type_filter(kwargs, elem_list)
        else:
            elems = iterator(kwargs, *args)

        super(Lattice, self).__init__(elems)

        # set default values
        kwargs.setdefault('name', '')
        kwargs.setdefault('periodicity', 1)
        if 'energy' not in kwargs:
            raise AtError('Lattice energy is not defined')
        # set attributes
        self.update(kwargs)

    def __getitem__(self, key):
        try:
            return super(Lattice, self).__getitem__(key.__index__())
        except (AttributeError, TypeError):
            return Lattice(self.iterator(key), iterator=no_filter, **vars(self))

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
        k1 = [(k, attrs.pop(k)) for k in Lattice._1st_attributes]
        k2 = [(k, v) for k, v in attrs.items() if not k.startswith('_')]
        keys = ', '.join('{0}={1!r}'.format(k, v) for k, v in (k1 + k2))
        return 'Lattice({0}, {1})'.format(super(Lattice, self).__repr__(), keys)

    def __str__(self):
        attrs = vars(self).copy()
        k1 = [(k, attrs.pop(k)) for k in Lattice._1st_attributes]
        k2 = [(k, v) for k, v in attrs.items() if not k.startswith('_')]
        keys = ', '.join('{0}={1!r}'.format(k, v) for k, v in (k1 + k2))
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

    def update(self, *args, **kwargs):
        """Update the element attributes with the given arguments

        update(**kwargs)
        update(mapping, **kwargs)
        update(iterable, **kwargs)
        """
        attrs = dict(*args, **kwargs)
        for (key, value) in attrs.items():
            setattr(self, key, value)

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
        try:
            return self._radiation
        except AttributeError:
            radiate = False
            for elem in self:
                if (elem.PassMethod.endswith('RadPass') or
                        elem.PassMethod.endswith('CavityPass')):
                    radiate = True
                    break
            # noinspection PyAttributeOutsideInit
            self._radiation = radiate
            return radiate

    # noinspection PyShadowingNames
    def modify_elements(self, elem_modify, copy=True):
        """Modify selected elements, in-place or in a lattice copy

        PARAMETERS
            elem_modify         element selection function.

            If elem_modify(elem) returns None, the element is unchanged.
            Otherwise, elem_modify(elem) must return a dictionary of
            attribute name and values, to be set to elem.

        RETURNS
            New lattice if copy == True
            None if copy == False

        KEYWORDS
            copy=True           Return a shallow copy of the lattice. Only the
                                modified attributes are replaced. Otherwise
                                modify the lattice in-place.
        """

        def lattice_modify():
            """Modifies the Lattice with elements modified by elem_modify"""
            radiate = False
            for elem in self:
                attrs = elem_modify(elem)
                if attrs is not None:
                    elem.update(attrs)
                if (elem.PassMethod.endswith('RadPass') or
                        elem.PassMethod.endswith('CavityPass')):
                    radiate = True
            self._radiation = radiate

        def lattice_copy(params):
            """Custom iterator for the creation of a new lattice"""
            radiate = False
            for elem in lattice_filter(params, self):
                attrs = elem_modify(elem)
                if attrs is not None:
                    elem = elem.copy()
                    elem.update(attrs)
                if (elem.PassMethod.endswith('RadPass') or
                        elem.PassMethod.endswith('CavityPass')):
                    radiate = True
                yield elem
            params['_radiation'] = radiate

        if copy:
            return Lattice(iterator=lattice_copy)
        else:
            lattice_modify()

    @staticmethod
    def _radiation_attrs(cavity_func, dipole_func,
                         quadrupole_func, wiggler_func,
                         sextupole_func, octupole_func):
        """Create a function returning the modified attributes"""

        def elem_func(elem):

            def isdipole(el):
                return isinstance(el, elements.Dipole) and (
                        el.BendingAngle != 0.0)

            if isinstance(elem, elements.RFCavity):
                return cavity_func(elem)
            elif isdipole(elem):
                return dipole_func(elem)
            elif isinstance(elem, elements.Quadrupole):
                return quadrupole_func(elem)
            elif isinstance(elem, elements.Wiggler):
                return wiggler_func(elem)
            elif isinstance(elem, elements.Sextupole):
                return sextupole_func(elem)
            elif isinstance(elem, elements.Octupole):
                return octupole_func(elem)
            else:
                return None

        return elem_func

    # noinspection PyShadowingNames
    def radiation_on(self, cavity_pass='CavityPass', dipole_pass='auto',
                     quadrupole_pass=None, wiggler_pass='auto',
                     sextupole_pass=None, octupole_pass=None, copy=False):
        """
        Turn acceleration and radiation on and return the lattice

        KEYWORDS
            cavity_pass='CavityPass'    PassMethod set on cavities
            dipole_pass='auto'          PassMethod set on dipoles
            quadrupole_pass=None        PassMethod set on quadrupoles
            wiggler_pass='auto'         PassMethod set on wigglers
            copy=False                  Return a shallow copy of the lattice and
                                        replace only the modified attributes
                                        Otherwise modify the lattice in-place.

            For PassMethod names, the convention is:
                None            no change
                'auto'          replace *Pass by *RadPass
                anything else   set as the new PassMethod
        """

        def repfunc(pass_method):
            if pass_method is None:
                # noinspection PyUnusedLocal
                def ff(elem):
                    return None
            elif pass_method == 'auto':
                def ff(elem):
                    if not elem.PassMethod.endswith('RadPass'):
                        pass_m = ''.join((elem.PassMethod[:-4], 'RadPass'))
                        return {'PassMethod': pass_m,
                                'Energy': self.energy}
                    else:
                        return None
            else:
                def ff(elem):
                    if elem.PassMethod != pass_method:
                        return {'PassMethod': pass_method,
                                'Energy': self.energy}
                    else:
                        return None
            return ff

        elem_func = self._radiation_attrs(repfunc(cavity_pass),
                                          repfunc(dipole_pass),
                                          repfunc(quadrupole_pass),
                                          repfunc(wiggler_pass),
                                          repfunc(sextupole_pass),
                                          repfunc(octupole_pass))
        return self.modify_elements(elem_func, copy=copy)

    # noinspection PyShadowingNames
    def radiation_off(self, cavity_pass='auto', dipole_pass='auto',
                      quadrupole_pass='auto', wiggler_pass='auto',
                      sextupole_pass='auto', octupole_pass='auto', copy=False):
        """
        Turn acceleration and radiation off and return the lattice

        KEYWORDS
            cavity_pass='IdentityPass'  PassMethod set on cavities
            dipole_pass='auto'          PassMethod set on dipoles
            quadrupole_pass=None        PassMethod set on quadrupoles
            wiggler_pass='auto'         PassMethod set on wigglers
            copy=False                  Return a shallow copy of the lattice and
                                        replace only the modified attributes
                                        Otherwise modify the lattice in-place.

            For PassMethod names, the convention is:
                None            no change
                'auto'          replace *RadPass by *Pass
                anything else   set as it is
        """

        def auto_cavity_pass(elem):
            newpass = 'IdentityPass' if elem.Length == 0 else 'DriftPass'
            if elem.PassMethod != newpass:
                return {'PassMethod': newpass}
            else:
                return None

        def auto_multipole_pass(elem):
            if elem.PassMethod.endswith('RadPass'):
                newpass = ''.join((elem.PassMethod[:-7], 'Pass'))
                return {'PassMethod': newpass}
            else:
                return None

        def repfunc(pass_method, auto_method):
            if pass_method is None:
                # noinspection PyUnusedLocal
                def ff(elem):
                    return None
            elif pass_method == 'auto':
                ff = auto_method
            else:
                def ff(elem):
                    if elem.PassMethod != pass_method:
                        return {'PassMethod': pass_method}
                    else:
                        return None
            return ff

        elem_func = self._radiation_attrs(
            repfunc(cavity_pass, auto_cavity_pass),
            repfunc(dipole_pass, auto_multipole_pass),
            repfunc(quadrupole_pass, auto_multipole_pass),
            repfunc(wiggler_pass, auto_multipole_pass),
            repfunc(sextupole_pass, auto_multipole_pass),
            repfunc(octupole_pass, auto_multipole_pass))
        return self.modify_elements(elem_func, copy=copy)

    def sbreak(self, break_s, break_elems=None):
        """Insert elements at selected locations in the lattice

        PARAMETERS
            break_s:        location or array of locations of breakpoints
            break_elems:    elements to be inserted at breakpoints (array of
                            elements as long as break_s or single element
                            duplicated as necessary). Default: Marker('sbreak')
        RETURNS
            A new lattice with new elements inserted at breakpoints
        """

        def sbreak_iterator(elems, insertions):
            """Iterate over elements and breaks where necessary"""

            def next_mk():
                """Extract the next element to insert"""
                try:
                    return next(insertions)
                except StopIteration:
                    return sys.float_info.max, None

            s_end = 0.0
            # get the 1st insertion
            smk, mk = next_mk()
            # skip all insertions at negative break_s, if any
            while smk < s_end:
                smk, mk = next_mk()

            for elem in elems:
                s_end += elem.Length
                # loop over all insertions within the element
                while smk < s_end:
                    frac = (s_end - smk) / elem.Length
                    if frac < 1.0:  # breakpoint is within the element
                        el0, elem = elem.divide([1.0 - frac, frac])
                        yield el0
                    yield mk
                    smk, mk = next_mk()
                yield elem

        # set default insertion
        if break_elems is None:
            break_elems = elements.Marker('sbreak')
        break_elems = numpy.reshape(break_elems, -1)
        # Check element lengths
        if not all(e.Length==0 for e in break_elems):
            warn(AtWarning(
                 "Inserting elements with length!=0 may change the lattice"))
        # broadcast break_s and break_elems to arrays of same size
        # and create an iterator over the elements to be inserted
        iter_mk = zip(*numpy.broadcast_arrays(break_s, break_elems))

        return Lattice(sbreak_iterator(self, iter_mk), **vars(self))


def lattice_filter(params, elems):
    """Copy lattice parameters an run through all lattice elements"""
    for key, value in vars(elems).items():
        params.setdefault(key, value)
    return iter(elems)


# noinspection PyUnusedLocal
def no_filter(params, elems):
    """Run through all elements without any check"""
    return iter(elems)


def type_filter(params, elems):
    """Run through all elements and check element validity.
    Analyses elements for radiation state
    """
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
    params['_radiation'] = radiate


def params_filter(params, elem_iterator, *args):
    """Run through all elements, looking for energy and periodicity.
    Remove the Energy attribute of non-radiating elements

    energy is taken from:
        1) The params dictionary
        2) Radiative elements (Cavity, RingParam, Radiating magnets)
        3) Any other element

    periodicity is taken from:
        1) The params dictionary
        2) Sum of the bending angles of magnets
    """
    rad_energies = []
    el_energies = []
    thetas = []

    for idx, elem in enumerate(elem_iterator(params, *args)):
        if (isinstance(elem, (elements.RFCavity, elements.Wiggler)) or
                elem.PassMethod.endswith('RadPass')):
            rad_energies.append(elem.Energy)
            try:
                elem.Energy = params['energy']
            except KeyError:
                params['energy'] = elem.Energy
        elif hasattr(elem, 'Energy'):
            el_energies.append(elem.Energy)
            del elem.Energy
        if isinstance(elem, elements.Dipole):
            thetas.append(elem.BendingAngle)
        yield elem

    energies = rad_energies or el_energies

    if 'energy' not in params and energies:
        # Guess energy from the Energy attribute of the elements
        params['energy'] = max(energies)

    if energies and min(energies) < max(energies):
        warn(AtWarning('Inconsistent energy values, '
                       '"energy" set to {0}'.format(params['energy'])))

    if 'periodicity' not in params:
        # Guess periodicity from the bending angles of the lattice
        try:
            nbp = 2.0 * numpy.pi / sum(thetas)
        except ZeroDivisionError:
            periodicity = 1
            warn(AtWarning('No bending in the cell, set "Periodicity" to 1'))
        else:
            periodicity = int(round(nbp))
            if abs(periodicity - nbp) > TWO_PI_ERROR:
                warn(AtWarning('Non-integer number of cells: '
                               '{0} -> {1}'.format(nbp, periodicity)))
        params['periodicity'] = periodicity


Lattice.get_elements = get_elements
Lattice.get_s_pos = get_s_pos
