"""Lattice object"""
import sys
import copy
import at  # for AtWarning, AtError
import numpy
from scipy.constants import physical_constants as cst
from warnings import warn
from math import pi
from . import elements
from .utils import checktype

__all__ = ['Lattice']

TWO_PI_ERROR = 1.E-4


class Lattice(list):
    """Lattice object
    An AT lattice is a sequence of AT elements. This object accepts extended indexing"""

    def __init__(self, *args, **kwargs):
        """Lattice(elements, **kwargs)
        Create a new lattice object

        INPUT
            elements:       iterable of AT elements

        KEYWORDS
            name:           name of the lattice (default '')
            energy:         energy of the lattice (default taken from the elements)
            periodicity:    numner of periods (default take from the elements)

            all other keywords will be set as Lattice attributes
        """
        descr = args[0] if len(args) > 0 else []
        if isinstance(descr, Lattice):
            params = []
            attrs = descr.__dict__
        else:
            params = [elem for elem in descr if isinstance(elem, elements.RingParam)]
            descr = [elem for elem in descr if not isinstance(elem, elements.RingParam)]
            radon = False
            for elem in descr:
                if elem.PassMethod.endswith('RadPass') or elem.PassMethod.endswith('CavityPass'):
                    radon = True
                    break
            attrs = dict(name='', energy=None, periodicity=None, _radiation_on=radon)
        attrs.update(kwargs)

        super(Lattice, self).__init__(descr)

        if attrs['energy'] is None:
            if len(params) > 0:
                attrs['energy'] = params[0].Energy
            else:
                # Guess energy from the Energy attribute of the elements
                energies = [elem.Energy for elem in descr if hasattr(elem, 'Energy')]
                if len(energies) == 0:
                    raise at.AtError('Energy not defined')
                energy = max(energies)
                if min(energies) < energy:
                    warn(at.AtWarning('Inconsistent energy values, "Energy" set to {0}'.format(energy)))
                attrs['energy'] = energy

        if attrs['periodicity'] is None:
            if len(params) > 0:
                attrs['periodicity'] = params[0].Periodicity
            else:
                # Guess periodicity from the bending angle of the superperiod
                theta = [elem.BendingAngle for elem in descr if isinstance(elem, elements.Dipole)]
                try:
                    nbp = 2.0 * pi / sum(theta)
                except ZeroDivisionError:
                    warn(at.AtWarning('No bending in the cell, set "Periodicity" to 1'))
                    attrs['periodicity'] = 1
                else:
                    periodicity = int(round(nbp))
                    if abs(periodicity - nbp) > TWO_PI_ERROR:
                        warn(at.AtWarning('non-integer number of cells: {0} -> {1}'.format(nbp, periodicity)))
                    attrs['periodicity'] = periodicity

        for key, value in attrs.items():
            setattr(self, key, value)

    def __getitem__(self, key):
        if isinstance(key, (int, numpy.int_)):
            return super(Lattice, self).__getitem__(key)
        elif isinstance(key, slice):
            return Lattice(super(Lattice, self).__getitem__(key), **self.__dict__)
        elif isinstance(key, numpy.ndarray) and key.dtype == bool:
            return Lattice([el for el, tst in zip(self, key) if tst], **self.__dict__)
        else:
            return Lattice([self[i] for i in key], **self.__dict__)

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
        volts = [elem.Voltage for elem in self if isinstance(elem, elements.RFCavity)]
        return self.periodicity * sum(volts)

    @property
    def harmonic_number(self):
        """Harmonic number"""
        harms = [elem.HarmNumber for elem in self if isinstance(elem, elements.RFCavity)]
        return self.periodicity * harms[0] if len(harms) > 0 else None

    @property
    def energy_loss(self):
        """Energy loss per turn

        Losses = Cgamma / 2pi * EGeV^4 * I2
        """
        lenthe = numpy.array([(elem.Length, elem.BendingAngle) for elem in self if isinstance(elem, elements.Dipole)])
        lendp = lenthe[:, 0]
        theta = lenthe[:, 1]

        e_radius = cst['classical electron radius'][0]
        e_mass = cst['electron mass energy equivalent in MeV'][0]
        cgamma = 4.0E9 * pi * e_radius / 3.0 / pow(e_mass, 3)

        I2 = self.periodicity * (numpy.sum(theta * theta / lendp))
        e_loss = cgamma / 2.0 / pi * pow(self.energy * 1.0E-9, 4) * I2 * 1.e9
        return e_loss

    def _radiation_switch(self, cavity_func=None, dipole_func=None, quadrupole_func=None):

        def mod_elem(ring, checkfun, modfun):
            n = 0
            for elem in filter(checkfun, ring):
                modfun(elem)
                n += 1
            return n

        if cavity_func is not None:
            n = mod_elem(self, checktype(elements.RFCavity), cavity_func)
            print('{0} modified cavities'.format(n))
        if dipole_func is not None:
            n = mod_elem(self, checktype(elements.Dipole), dipole_func)
            print('{0} modified dipoles'.format(n))
        if quadrupole_func is not None:
            n = mod_elem(self, checktype(elements.Quadrupole), quadrupole_func)
            print('{0} modified quadrupoles'.format(n))

    def radiation_on(self, cavity_pass='CavityPass', dipole_pass='auto', quadrupole_pass=None):
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
                        elem.PassMethod = ''.join((elem.PassMethod[:-4], 'RadPass'))
                        elem.Energy = self.energy
            else:
                def ff(elem):
                    elem.PassMethod = pass_method
                    elem.Energy = self.energy
            return ff

        self._radiation_switch(repfunc(cavity_pass), repfunc(dipole_pass), repfunc(quadrupole_pass))
        self._radiation_on = True

    def radiation_off(self, cavity_pass='IdentityPass', dipole_pass='auto', quadrupole_pass=None):
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
                        elem.PassMethod = ''.join((elem.PassMethod[:-7], 'Pass'))
            else:
                def ff(elem):
                    elem.PassMethod = pass_method
            return ff

        self._radiation_switch(repfunc(cavity_pass), repfunc(dipole_pass), repfunc(quadrupole_pass))
        self._radiation_on = False
