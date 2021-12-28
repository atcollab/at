"""Lattice object

The methods implemented in this module are internal to the 'lattice' package.
This is necessary to ensure that the 'lattice' package is independent of other
AT packages.

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
from .constants import clight
from .particle_object import Particle
from .utils import AtError, AtWarning
from .utils import uint32_refpts as uint32_refs, bool_refpts as bool_refs
from .utils import refpts_iterator, refpts_len
from .utils import get_s_pos, get_elements, get_cells, get_refpts
from .utils import get_value_refpts, set_value_refpts
from .utils import set_shift, set_tilt
from . import elements
# noinspection PyProtectedMember
from .utils import _uint32_refs, _bool_refs

TWO_PI_ERROR = 1.E-4

__all__ = ['Lattice', 'type_filter', 'params_filter', 'lattice_filter',
           'no_filter']

# Don't warn on floating-pont errors
numpy.seterr(divide='ignore', invalid='ignore')


# noinspection PyAttributeOutsideInit
class Lattice(list):
    """Lattice object
    An AT lattice is a sequence of AT elements.
    A Lattice accepts extended indexing (as a numpy ndarray).

    Lattice attributes;
        name            Name of the lattice
        energy          Particle energy
        periodicity     Number of super-periods to describe the full ring
        particle        Circulating particle
        harmonic_number Harmonic number of the full ring (periodicity x cells)

    Lattice(elems, **params)        Create a new lattice object

        INPUT
            elems:                  any iterable of AT elements

        KEYWORDS
            name=''                 Name of the lattice
            energy                  Energy of the lattice
            periodicity=1           Number of periods
            particle='relativistic' Circulating particle. May be 'relativistic',
                                    'electron', 'positron', 'proton'
                                    or a Particle object
            iterator=None           Custom iterator (see below)
            *                       All other keywords will be set as attributes
                                    of the Lattice object

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
                   with the keywords of the lattice constructor. The custom
                   iterator may add, remove or mofify parameters.
                   Finally, the remaining parameters will be set as Lattice
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
    # Attributes displayed:
    _disp_attributes = ('name', 'energy', 'particle', 'periodicity',
                        'harmonic_number')
    # Attributes propagated in copies:
    _std_attributes = ('name', '_energy', '_particle', 'periodicity',
                       '_cell_harmnumber', '_radiation')

    # noinspection PyUnusedLocal
    def __init__(self, *args, iterator=None, scan=False, **kwargs):
        """Lattice constructor"""

        if iterator is None:
            arg1, = args or [[]]  # accept 0 or 1 argument
            if isinstance(arg1, Lattice):
                elems = arg1.attrs_filter(kwargs, arg1)
            else:
                elems = params_filter(kwargs, type_filter, arg1)
        else:
            elems = iterator(kwargs, *args)

        super(Lattice, self).__init__(elems)

        # set default values
        kwargs.setdefault('name', '')
        periodicity = kwargs.setdefault('periodicity', 1)
        kwargs.setdefault('_particle', Particle('electron'))
        # Remove temporary keywords
        frequency = kwargs.pop('_frequency', None)
        cell_h = kwargs.pop('_harmnumber', None)

        if 'harmonic_number' in kwargs:
            cell_h = kwargs.pop('harmonic_number') / periodicity
        if 'energy' in kwargs:
            kwargs.pop('_energy', None)
        elif '_energy' not in kwargs:
            raise AtError('Lattice energy is not defined')
        if 'particle' in kwargs:
            kwargs.pop('_particle', None)
        # set attributes
        self.update(kwargs)

        if cell_h is not None:
            self._cell_harmnumber = cell_h
        elif frequency is not None:
            gamma = self.gamma
            beta = math.sqrt(1.0 - 1.0 / gamma / gamma)
            rev = beta * clight / self.get_s_pos(len(self))[0]
            self._cell_harmnumber = int(round(frequency / rev))

    def __getitem__(self, key):
        try:                                # Integer
            return super(Lattice, self).__getitem__(key.__index__())
        except (AttributeError, TypeError):
            if isinstance(key, slice):      # Slice
                rg = range(*key.indices(len(self)))
            else:                           # Array of integers or boolean
                rg = self.uint32_refpts(key)
            return Lattice((super(Lattice, self).__getitem__(i) for i in rg),
                           iterator=self.attrs_filter)

    def __setitem__(self, key, values):
        try:                                # Integer or slice
            super(Lattice, self).__setitem__(key, values)
        except TypeError:                   # Array of integers or boolean
            rg = self.uint32_refpts(key)
            for i, v in zip(*numpy.broadcast_arrays(rg, values)):
                super(Lattice, self).__setitem__(i, v)

    def __delitem__(self, key):
        try:                                # Integer or slice
            super(Lattice, self).__delitem__(key)
        except TypeError:                   # Array of integers or boolean
            rg = self.uint32_refpts(key)
            for i in reversed(rg):
                super(Lattice, self).__delitem__(i)

    def __repr__(self):
        ks = ', '.join('{0}={1!r}'.format(k, v) for k, v in self.attrs.items())
        return 'Lattice({0}, {1})'.format(super(Lattice, self).__repr__(), ks)

    def __str__(self):
        ks = ', '.join('{0}={1!r}'.format(k, v) for k, v in self.attrs.items())
        return 'Lattice(<{0} elements>, {1})'.format(len(self), ks)

    def __add__(self, elems):
        """Add elems, an iterable of AT elements, to the lattice"""
        def add_filter(params, el1, el2):
            it1 = el1.attrs_filter(params, el1)
            return type_filter(params, itertools.chain(it1, el2))
        return Lattice(self, elems, iterator=add_filter)

    def __mul__(self, n):
        """Repeats n times the lattice"""
        return Lattice(itertools.chain(*itertools.repeat(self, n)),
                       iterator=self.attrs_filter)

    @property
    def attrs(self):
        """Dictionary of lattice attributes"""
        def extattr(d):
            for k in self._disp_attributes:
                d.pop(k, None)
                yield k, getattr(self, k, None)

        vrs = vars(self).copy()
        # Standard attributes
        res = {k: v for k, v in extattr(vrs) if v is not None}
        # Custom attributes
        res.update((k, v) for k, v in vrs.items() if not k.startswith('_'))
        return res

    def uint32_refpts(self, refpts):
        """"Return a uint32 numpy array containing the indices of the selected
        elements
        """
        if callable(refpts):
            refpts = [refpts(el) for el in self]
        return uint32_refs(refpts, len(self))

    def bool_refpts(self, refpts):
        """Return a boolean numpy array of length n_elements + 1 where
        True elements are selected.
        """
        if callable(refpts):
            refpts = [refpts(el) for el in self]
        return bool_refs(refpts, len(self))

    def rotate(self, n):
        """Return a new lattice rotated left by n elements"""
        if len(self) == 0:
            return self.copy()
        n = n % len(self)      # works for n<0
        return self[n:] + self[:n]

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
        """Return a shallow copy of the lattice"""
        return copy.copy(self)

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
            yield from self[:ibeg]
            for elem in self[ibeg:iend]:
                nslices = int(math.ceil(elem.Length / size))
                if nslices > 1:
                    frac = numpy.ones(nslices) / nslices
                    for el in elem.divide(frac):
                        yield el
                else:
                    yield elem
            yield from self[iend:]

        s_range = self.s_range
        if size is None:
            smin, smax = s_range
            size = (smax - smin) / slices
        i_range = self.i_range
        return Lattice(slice_iter(i_range[0], i_range[-1]),
                       iterator=self.attrs_filter, s_range=s_range)

    @property
    def attrs_filter(self):
        """Filter function which duplicates the lattice attributes"""
        def filt(params, elems_iterator):
            for key in self._std_attributes:
                try:
                    params.setdefault(key, getattr(self, key))
                except AttributeError:
                    pass
            return elems_iterator

        return filt

    @property
    def s_range(self):
        """Range of interest: [s_min, s_max]. 'None' means the full cell."""
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
        return self.uint32_refpts(i_range)

    @property
    def energy(self):
        """Lattice energy"""
        return self._energy

    @energy.setter
    def energy(self, energy):
        # Set the Energy attribute of radiating elements
        for elem in self:
            if (isinstance(elem, (elements.RFCavity, elements.Wiggler)) or
                    elem.PassMethod.endswith('RadPass')):
                elem.Energy = energy
        # Set the energy attribute of the Lattice
        # Use a numpy scalar to allow division by zero
        # self._energy = numpy.array(energy, dtype=float)
        self._energy = energy

    @property
    def circumference(self):
        """Ring circumference (full ring) [m]"""
        return self.periodicity * self.get_s_pos(len(self))[0]

    @property
    def revolution_frequency(self):
        """Revolution frequency (fullring) [Hz]"""
        # gamma = self.gamma
        # beta = math.sqrt(1.0 - 1.0 / gamma / gamma)
        # return beta * clight / self.circumference
        return clight / self.circumference

    @property
    def particle(self):
        """Circulating particle"""
        return self._particle

    @particle.setter
    def particle(self, particle):
        if isinstance(particle, Particle):
            self._particle = particle
        else:
            self._particle = Particle(particle)

    @property
    def harmonic_number(self):
        try:
            return self.periodicity * self._cell_harmnumber
        except AttributeError:
            raise AttributeError('harmonic_number undefined: '
                                 'No cavity found in the lattice') from None

    @harmonic_number.setter
    def harmonic_number(self, value):
        periods = int(self.periodicity)
        cell_h, rem = divmod(int(value), periods)
        if rem == 0:
            self._cell_harmnumber = cell_h
        else:
            raise AtError('harmonic number must be a multiple of {}'
                          .format(periods))

    @property
    def gamma(self):
        return float(self.energy / self.particle.rest_energy)

    @property
    def beta(self):
        gamma = float(self.energy / self.particle.rest_energy)
        return math.sqrt(1.0 - 1.0/gamma/gamma)

    # noinspection PyPep8Naming
    @property
    def BRho(self):
        return math.sqrt(self.energy**2 - self.particle.rest_energy**2)/clight

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
    def modify_elements(self, elem_modify, copy=True, **kwargs):
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
            copy=True   If True, return a shallow copy of the lattice. Only the
                        modified elements are copied.
                        If False, the modification is done in-place
        """

        def lattice_modify(**kws):
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
            self.update(kws)

        def lattice_copy(params):
            """Custom iterator for the creation of a new lattice"""
            radiate = False
            for elem in self.attrs_filter(params, self):
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
            return Lattice(iterator=lattice_copy, **kwargs)
        else:
            lattice_modify(**kwargs)

    @staticmethod
    def _radiation_attrs(cavity_func, dipole_func,
                         quadrupole_func, wiggler_func,
                         sextupole_func, octupole_func, 
                         multipole_func):
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
            elif isinstance(elem, elements.Multipole):
                return multipole_func(elem)
            else:
                return None

        return elem_func

    # noinspection PyShadowingNames
    def radiation_on(self, cavity_pass='CavityPass', dipole_pass='auto',
                     quadrupole_pass='auto', wiggler_pass='auto',
                     sextupole_pass=None, octupole_pass=None, 
                     multipole_pass=None, copy=False):
        """
        Turn acceleration and radiation on and return the lattice

        KEYWORDS
            cavity_pass='CavityPass'    PassMethod set on cavities
            dipole_pass='auto'          PassMethod set on dipoles
            quadrupole_pass='auto'      PassMethod set on quadrupoles
            wiggler_pass='auto'         PassMethod set on wigglers
            copy=False  If False, the modification is done in-place,
                        If True, return a shallow copy of the lattice. Only the
                        radiating elements are copied with PassMethod modified.
                        CAUTION: a shallow copy means that all non-radiating
                        elements are shared with the original lattice.
                        Any further modification will affect in both lattices.

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
                                          repfunc(octupole_pass),
                                          repfunc(multipole_pass))
        return self.modify_elements(elem_func, copy=copy)

    # noinspection PyShadowingNames
    def radiation_off(self, cavity_pass='auto', dipole_pass='auto',
                      quadrupole_pass='auto', wiggler_pass='auto',
                      sextupole_pass='auto', octupole_pass='auto', 
                      multipole_pass='auto', copy=False):
        """
        Turn acceleration and radiation off and return the lattice

        KEYWORDS
            cavity_pass='IdentityPass'  PassMethod set on cavities
            dipole_pass='auto'          PassMethod set on dipoles
            quadrupole_pass=None        PassMethod set on quadrupoles
            wiggler_pass='auto'         PassMethod set on wigglers
            copy=False  If False, the modification is done in-place,
                        If True, return a shallow copy of the lattice. Only the
                        radiating elements are copied with PassMethod modified.
                        CAUTION: a shallow copy means that all non-radiating
                        elements are shared with the original lattice.
                        Any further modification will affect in both lattices.

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
            repfunc(octupole_pass, auto_multipole_pass),
            repfunc(multipole_pass, auto_multipole_pass))
        return self.modify_elements(elem_func, copy=copy)

    def sbreak(self, break_s, break_elems=None, **kwargs):
        """Insert elements at selected locations in the lattice

        PARAMETERS
            break_s:        location or array of locations of breakpoints
            break_elems:    elements to be inserted at breakpoints (array of
                            elements as long as break_s or single element
                            duplicated as necessary). Default: Marker('sbreak')
        RETURNS
            A new lattice with new elements inserted at breakpoints
        """

        def sbreak_iterator(params):
            """Iterate over elements and breaks where necessary"""

            def next_mk():
                """Extract the next element to insert"""
                try:
                    return next(iter_mk)
                except StopIteration:
                    return sys.float_info.max, None

            s_end = 0.0
            # get the 1st insertion
            smk, mk = next_mk()
            # skip all insertions at negative break_s, if any
            while smk < s_end:
                smk, mk = next_mk()

            for elem in self.attrs_filter(params, self):
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
        if not all(e.Length == 0 for e in break_elems):
            warn(AtWarning(
                 "Inserting elements with length!=0 may change the lattice"))
        # broadcast break_s and break_elems to arrays of same size
        # and create an iterator over the elements to be inserted
        iter_mk = zip(*numpy.broadcast_arrays(break_s, break_elems))

        return Lattice(iterator=sbreak_iterator, **kwargs)

    def replace(self, refpts, **kwargs):
        """Return a shallow copy of the lattice replacing the selected
        elements by a deep copy"""
        if callable(refpts):
            check = map(refpts, self)
        else:
            check = iter(self.bool_refpts(refpts))
        elems = (el.deepcopy() if ok else el for el, ok in zip(self, check))
        return Lattice(elems, iterator=self.attrs_filter, **kwargs)


def lattice_filter(params, lattice):
    """Copy lattice parameters and run through all lattice elements"""
    return lattice.attrs_filter(params, lattice)


# noinspection PyUnusedLocal
def no_filter(params, elems):
    """Run through all elements without any check"""
    return elems


def type_filter(params, elem_iterator):
    """Run through all elements and check element validity.
    Analyse elements for radiation state
    """
    radiate = False
    for idx, elem in enumerate(elem_iterator):
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
        2) Cavity elements
        3) Any other element

    periodicity is taken from:
        1) The params dictionary
        2) Sum of the bending angles of magnets
    """
    el_energies = []
    thetas = []
    cavities = []

    for idx, elem in enumerate(elem_iterator(params, *args)):
        if isinstance(elem, elements.RFCavity):
            cavities.append(elem)
        elif hasattr(elem, 'Energy'):
            el_energies.append(elem.Energy)
            del elem.Energy
        if isinstance(elem, elements.Dipole):
            thetas.append(elem.BendingAngle)
        yield elem

    cav_energies = [el.Energy for el in cavities if hasattr(el, 'Energy')]
    cavities.sort(key=lambda el: el.Frequency)
    if cavities:
        params['_harmnumber'] = getattr(cavities[0], 'HarmNumber', None)
        params['_frequency'] = getattr(cavities[0], 'Frequency', None)

    if 'energy' not in params:
        energies = cav_energies or el_energies
        if energies:
            # Guess energy from the Energy attribute of the elements
            energy = params.setdefault('energy',  max(energies))
            if min(energies) < max(energies):
                warn(AtWarning('Inconsistent energy values, '
                               '"energy" set to {0}'.format(energy)))

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


Lattice.uint32_refpts = _uint32_refs
Lattice.bool_refpts = _bool_refs
Lattice.get_cells = get_cells
Lattice.get_refpts = get_refpts
Lattice.set_shift = set_shift
Lattice.set_tilt = set_tilt
Lattice.get_elements = get_elements
Lattice.get_s_pos = get_s_pos
Lattice.select = refpts_iterator
Lattice.refcount = refpts_len
Lattice.get_value_refpts = get_value_refpts
Lattice.set_value_refpts = set_value_refpts
