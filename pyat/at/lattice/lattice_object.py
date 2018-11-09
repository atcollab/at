"""Lattice object"""
import sys
import copy
import numpy
from scipy.constants import physical_constants as cst
from warnings import warn
from math import pi
from . import elements, checktype, AtWarning, AtError
from ..lattice import get_s_pos
from ..physics import find_orbit4, find_orbit6, find_sync_orbit, find_m44, \
    find_m66, linopt, ohmi_envelope, get_mcf

__all__ = ['Lattice']

TWO_PI_ERROR = 1.E-4


def _scan_list(descr, keep_all=False, **kwargs):
    """Extract the lattice attributes from an element list"""
    attrs = {'name': '', '_radiation_on': False}
    params = [elem for elem in descr if isinstance(elem, elements.RingParam)]
    if params:
        # Set all RingParam attributes to the lattice object
        attrs.update(
            (Lattice._translate.get(key, key.lower()), value)
            for (key, value) in vars(params[0]).items()
            if key not in Lattice._ignore)

        if not keep_all:
            descr = [elem for elem in descr
                     if not isinstance(elem, elements.RingParam)]

    for elem in descr:
        if elem.PassMethod.endswith('RadPass') or elem.PassMethod.endswith(
                'CavityPass'):
            attrs['_radiation_on'] = True
            break

    attrs.update(kwargs)

    if attrs.get('energy') is None:
        # Guess energy from the Energy attribute of the elements
        # Look first in cavities
        energies = [elem.Energy for elem in descr
                    if isinstance(elem, elements.RFCavity)]
        if not energies:
            # Then look in all elements
            energies = [elem.Energy for elem in descr
                        if hasattr(elem, 'Energy')]
        if not energies:
            raise AtError('Lattice Energy is not defined')
        energy = max(energies)
        if min(energies) < energy:
            warn(AtWarning(
                'Inconsistent energy values, "Energy" set to {0}'.format(
                    energy)))
        attrs['energy'] = energy

    if attrs.get('periodicity') is None:
        # Guess periodicity from the bending angle of the superperiod
        theta = [elem.BendingAngle for elem in descr if
                 isinstance(elem, elements.Dipole)]
        try:
            nbp = 2.0 * pi / sum(theta)
        except ZeroDivisionError:
            warn(AtWarning('No bending in the cell, set "Periodicity" to 1'))
            attrs['periodicity'] = 1
        else:
            periodicity = int(round(nbp))
            if abs(periodicity - nbp) > TWO_PI_ERROR:
                warn(AtWarning(
                    'Non-integer number of cells: {0} -> {1}'.format(
                        nbp, periodicity)))
            attrs['periodicity'] = periodicity

    return descr, attrs


class Lattice(list):
    """Lattice object
    An AT lattice is a sequence of AT elements.
    This object accepts extended indexing"""

    _translate = {'Energy': 'energy', 'Periodicity': 'periodicity',
                  'FamName': 'name'}
    _ignore = {'PassMethod', 'Length'}

    def __init__(self, descr=None, _scan=True, **kwargs):
        """Lattice(elements, **kwargs)
        Create a new lattice object

        INPUT
            elements:       sequence of AT elements

        KEYWORDS
            keep_all=False  if True, keep RingParam elements in the lattice
            name=''         Name of the lattice
            energy          Energy of the lattice
                            (default: taken from the elements)
            periodicity     Number of periods
                            (default: taken from the elements)

            all other keywords will be set as Lattice attributes
        """
        if descr is None:
            descr = []

        if isinstance(descr, Lattice):
            attrs = {}
            attrs.update(vars(descr))
            attrs.update(kwargs)
        elif _scan:
            descr, attrs = _scan_list(descr, **kwargs)
        else:
            attrs = {}
            attrs.update(kwargs)

        super(Lattice, self).__init__(descr)

        for key, value in attrs.items():
            setattr(self, key, value)

    def __getitem__(self, key):
        try:
            elems = super(Lattice, self).__getitem__(key)
        except TypeError:
            if isinstance(key, numpy.ndarray) and key.dtype == bool:
                return Lattice((el for el, tst in zip(self, key) if tst),
                               _scan=False, **vars(self))
            else:
                return Lattice((self[i] for i in key), _scan=False,
                               **vars(self))
        else:
            if isinstance(elems, list):
                elems = Lattice(elems, _scan=False, **vars(self))
            return elems

    def __repr__(self):
        at = ', '.join(
            '{0}={1!r}'.format(key, val) for key, val in vars(self).items() if
            not key.startswith('_'))
        return 'Lattice({0}, {1})'.format(super(Lattice, self).__repr__(), at)

    def __str__(self):
        at = ', '.join(
            '{0}={1!r}'.format(key, val) for key, val in vars(self).items() if
            not key.startswith('_'))
        return 'Lattice(<{0} elements>, {1})'.format(len(self), at)

    def __add__(self, other):
        return Lattice(super(Lattice, self).__add__(other), _scan=False,
                       **vars(self))

    def __mul__(self, other):
        return Lattice(super(Lattice, self).__mul__(other), _scan=False,
                       **vars(self))

    if sys.version_info < (3, 0):
        # This won't be defined if version is at least 3.0
        def __getslice__(self, i, j):
            return self[max(0, i):max(0, j):]

    def copy(self):
        """Return a shallow copy"""
        return Lattice(self)

    def deepcopy(self):
        """Return a deep copy"""
        return copy.deepcopy(self)

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
        cgamma = 4.0E9 * pi * e_radius / 3.0 / pow(e_mass, 3)

        i2 = self.periodicity * (numpy.sum(theta * theta / lendp))
        e_loss = cgamma / 2.0 / pi * pow(self.energy * 1.0E-9, 4) * i2 * 1.e9
        return e_loss

    def _radiation_switch(self, cavity_func=None, dipole_func=None,
                          quadrupole_func=None):

        def mod_elem(ring, checkfun, modfun):
            nb = 0
            for elem in filter(checkfun, ring):
                modfun(elem)
                nb += 1
            return nb

        def checkdipole(elem):
            return isinstance(elem, elements.Dipole) and (
                    elem.BendingAngle != 0.0)

        if cavity_func is not None:
            n = mod_elem(self, checktype(elements.RFCavity), cavity_func)
            print('{0} modified cavities'.format(n))
        if dipole_func is not None:
            n = mod_elem(self, checkdipole, dipole_func)
            print('{0} modified dipoles'.format(n))
        if quadrupole_func is not None:
            n = mod_elem(self, checktype(elements.Quadrupole), quadrupole_func)
            print('{0} modified quadrupoles'.format(n))

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
        self._radiation_on = True

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
        self._radiation_on = False

    def linopt(self, *args, **kwargs):
        """See at.physics.linopt():
        """
        if self._radiation_on:
            raise AtError('linopt needs no radiation in the lattice')
        return linopt(self, *args, **kwargs)

    def ohmi_envelope(self, *args, **kwargs):
        """See at.physics.ohmi_envelope():
        """
        if not self._radiation_on:
            raise AtError('ohmi_envelope needs radiation in the lattice')
        return ohmi_envelope(self, *args, **kwargs)


Lattice.get_s_pos = get_s_pos
Lattice.get_mcf = get_mcf
Lattice.find_orbit4 = find_orbit4
Lattice.find_sync_orbit = find_sync_orbit
Lattice.find_orbit6 = find_orbit6
Lattice.find_m44 = find_m44
Lattice.find_m66 = find_m66

if sys.version_info < (3, 0):
    Lattice.linopt.__func__.__doc__ += linopt.__doc__
    Lattice.ohmi_envelope.__func__.__doc__ += ohmi_envelope.__doc__
else:
    Lattice.linopt.__doc__ += linopt.__doc__
    Lattice.ohmi_envelope.__doc__ += ohmi_envelope.__doc__
