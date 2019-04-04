"""Lattice object"""
import sys
import copy
import numpy
import math
from scipy.constants import physical_constants as cst
from warnings import warn
from at.lattice import elements, get_s_pos, get_elements, uint32_refpts
from at.lattice import AtWarning, AtError
from at.physics import find_orbit4, find_orbit6, find_sync_orbit, find_m44
from at.physics import find_m66, linopt, ohmi_envelope, get_mcf
from at.plot import plot_beta, plot_trajectory

__all__ = ['Lattice', 'get_ring_properties']

TWO_PI_ERROR = 1.E-4


class Lattice(list):
    """Lattice object
    An AT lattice is a sequence of AT elements.
    A Lattice accepts extended indexing (as a numpy ndarray).

    ATTRIBUTES
        name        Name of the lattice
        energy      Particle energy
        periodicity Number of super-periods to describe the full ring

        """

    _translate = {'Energy': 'energy', 'Periodicity': 'periodicity',
                  'FamName': 'name'}
    _ignore = {'PassMethod', 'Length'}

    @staticmethod
    def _scanner(elems, attributes, keep_all=False, **kwargs):
        """Extract the lattice attributes from an element list

        If a RingParam element is found, its energy will be used, energies
        defined in other elements are ignored. All its user-defined
        attributes will also be set as Lattice attributes.

        Otherwise, necessary values will be guessed from other elements.
        """
        params = []
        rf_energies = []
        el_energies = []
        thetas = []

        radiate = False
        for idx, elem in enumerate(elems):
            if isinstance(elem, elements.RingParam):
                params.append(elem)
                if keep_all:
                    yield elem
            elif isinstance(elem, elements.Element):
                yield elem
            else:
                warn(AtWarning('item {0} ({1}) is not an AT element: '
                               'ignored'.format(idx, elem)))
                continue

            if hasattr(elem, 'Energy'):
                if isinstance(elem, elements.RFCavity):
                    rf_energies.append(elem.Energy)
                el_energies.append(elem.Energy)
            if isinstance(elem, elements.Dipole):
                thetas.append(elem.BendingAngle)
            if elem.PassMethod.endswith(
                    'RadPass') or elem.PassMethod.endswith('CavityPass'):
                radiate = True
        attributes['_radiation'] = radiate

        if params:
            # At least one RingParam element, use the 1st one
            if len(params) > 1:
                warn(AtWarning(
                    'More than 1 RingParam element, the 1st one is used'))
            attributes.update((Lattice._translate.get(key, key.lower()),
                               getattr(params[0], key))
                              for key in params[0]
                              if key not in Lattice._ignore)
        else:
            # No RingParam element, try to guess
            attributes['name'] = ''
            if 'energy' not in kwargs:
                # Guess energy from the Energy attribute of the elements
                energies = rf_energies if rf_energies else el_energies
                if not energies:
                    raise AtError('Lattice energy is not defined')
                energy = max(energies)
                if min(energies) < energy:
                    warn(AtWarning('Inconsistent energy values, '
                                   '"energy" set to {0}'.format(energy)))
                attributes['energy'] = energy
            if 'periodicity' not in kwargs:
                # Guess periodicity from the bending angles of the lattice
                try:
                    nbp = 2.0 * numpy.pi / sum(thetas)
                except ZeroDivisionError:
                    warn(AtWarning('No bending in the cell, '
                                   'set "Periodicity" to 1'))
                    attributes['periodicity'] = 1
                else:
                    periodicity = int(round(nbp))
                    if abs(periodicity - nbp) > TWO_PI_ERROR:
                        warn(AtWarning(
                            'Non-integer number of cells: '
                            '{0} -> {1}'.format(nbp, periodicity)))
                    attributes['periodicity'] = periodicity

    def __init__(self, elems=None, name=None, energy=None, periodicity=None,
                 keep_all=False, **kwargs):
        """Create a new lattice object

        INPUT
            elems:          any iterable of AT elements

        KEYWORDS
            name            Name of the lattice
                            (default: taken from the lattice, or '')
            energy          Energy of the lattice
                            (default: taken from the elements)
            periodicity     Number of periods
                            (default: taken from the elements, or 1)
            keep_all=False  if True, keep RingParam elements in the lattice

            all other keywords will be set as Lattice attributes
        """
        if elems is None:
            elems = []

        if name is not None:
            kwargs['name'] = name
        if energy is not None:
            kwargs['energy'] = energy
        if periodicity is not None:
            kwargs['periodicity'] = periodicity

        attrs = {}
        if isinstance(elems, Lattice):
            attrs.update(vars(elems))
        elif None in {name, energy, periodicity}:
            elems = Lattice._scanner(iter(elems), attrs, keep_all=keep_all,
                                     **kwargs)

        super(Lattice, self).__init__(elems)

        attrs.update(kwargs)

        if '_radiation' not in attrs:
            attrs['_radiation'] = False
            for elem in self:
                if (elem.PassMethod.endswith('RadPass') or
                        elem.PassMethod.endswith('CavityPass')):
                    attrs['_radiation'] = True
                    break

        for key, value in attrs.items():
            setattr(self, key, value)

        self.s_range = getattr(self, '_s_range', None)

    def __getitem__(self, key):
        try:
            return super(Lattice, self).__getitem__(key)
        except TypeError:
            key = uint32_refpts(key, len(self))
            return [super(Lattice, self).__getitem__(i) for i in key]

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
        keywords = ['{0}={1!r}'.format(key, attrs.pop(key)) for key in
                    Lattice._translate.values()]
        keywords += ['{0}={1!r}'.format(key, val) for key, val in
                     attrs.items() if not key.startswith('_')]
        return 'Lattice({0}, {1})'.format(super(Lattice, self).__repr__(),
                                          ', '.join(keywords))

    def __str__(self):
        attrs = vars(self).copy()
        keywords = ['{0}={1!r}'.format(key, attrs.pop(key)) for key in
                    Lattice._translate.values()]
        keywords += ['{0}={1!r}'.format(key, val) for key, val in
                     attrs.items() if not key.startswith('_')]
        return 'Lattice(<{0} elements>, {1})'.format(len(self),
                                                     ', '.join(keywords))

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

        if size is None:
            smin, smax = self.s_range
            size = (smax - smin) / slices

        i1 = self._i_range[0]
        i2 = self._i_range[-1]
        return Lattice(slice_iter(i1, i2), **vars(self))

    @property
    def s_range(self):
        """Range of interest. Initially set to the full cell.
        'None' means the full cell."""
        return self._s_range

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
        return uint32_refpts(self._i_range, len(self))

    @property
    def voltage(self):
        """Accelerating voltage"""
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

    @property
    def energy_loss(self):
        """Energy loss per turn [eV]

        Losses = Cgamma / 2pi * EGeV^4 * i2
        """
        lenthe = numpy.array(
            [(elem.Length, elem.BendingAngle) for elem in self if
             isinstance(elem, elements.Dipole)])
        lendp = lenthe[:, 0]
        theta = lenthe[:, 1]

        e_radius = cst['classical electron radius'][0]
        e_mass = cst['electron mass energy equivalent in MeV'][0]
        cgamma = 4.0E9 * numpy.pi * e_radius / 3.0 / pow(e_mass, 3)

        i2 = self.periodicity * (numpy.sum(theta * theta / lendp))
        e_loss = cgamma / 2.0 / numpy.pi * pow(self.energy * 1.0E-9,
                                               4) * i2 * 1.e9
        return e_loss

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

    def linopt(self, *args, **kwargs):
        """See at.physics.linopt():
        """
        if self._radiation:
            raise AtError('linopt needs no radiation in the lattice')
        return linopt(self, *args, **kwargs)

    def get_mcf(self, *args, **kwargs):
        """See at.physics.get_mcf():
        """
        if self._radiation:
            raise AtError('get_mcf needs no radiation in the lattice')
        return get_mcf(self, *args, **kwargs)

    def ohmi_envelope(self, *args, **kwargs):
        """See at.physics.ohmi_envelope():
        """
        if not self._radiation:
            raise AtError('ohmi_envelope needs radiation in the lattice')
        return ohmi_envelope(self, *args, **kwargs)


def get_ring_properties(ring):
    """Get ring properties from a Lattice object or from a sequence of elements.

    PARAMETER
        ring            Lattice object or any iterable of AT elements

    RETURN
        Directory with keys:

            name        Name of the lattice taken, in priority order, from:
                            - a Lattice object,
                            - a RingParam element
                        Defaults to ''
            energy      Particle energy taken, in priority order, from:
                            - a Lattice object,
                            - a RingParam element
                            - a RFCavity element
                            - any other element having an 'Energy' attribute
                        An error is raised if no energy is found
            periodicity Number of super-periods to describe the full ring,
                        taken in priority order from:
                            - a Lattice object,
                            - a RingParam element
                            - the sum of bending angles in the lattice
                        Defaults to 1

        Other properties may come from a lattice object or a RingParam element
        in the lattice
    """
    if isinstance(ring, Lattice):
        props = vars(ring).copy()
    else:
        props = {}
        # noinspection PyUnusedLocal,PyProtectedMember
        for elem in Lattice._scanner(ring, props):
            pass
    return props


Lattice.get_elements = get_elements
Lattice.get_s_pos = get_s_pos
Lattice.find_orbit4 = find_orbit4
Lattice.find_sync_orbit = find_sync_orbit
Lattice.find_orbit6 = find_orbit6
Lattice.find_m44 = find_m44
Lattice.find_m66 = find_m66
Lattice.plot_beta = plot_beta
Lattice.plot_trajectory = plot_trajectory

if sys.version_info < (3, 0):
    # noinspection PyUnresolvedReferences
    Lattice.linopt.__func__.__doc__ += linopt.__doc__
    # noinspection PyUnresolvedReferences
    Lattice.get_mcf.__func__.__doc__ += get_mcf.__doc__
    # noinspection PyUnresolvedReferences
    Lattice.ohmi_envelope.__func__.__doc__ += ohmi_envelope.__doc__
else:
    Lattice.linopt.__doc__ += linopt.__doc__
    Lattice.get_mcf.__doc__ += get_mcf.__doc__
    Lattice.ohmi_envelope.__doc__ += ohmi_envelope.__doc__
