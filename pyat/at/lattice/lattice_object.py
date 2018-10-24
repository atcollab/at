"""Lattice object"""
import sys
import copy
import numpy
from scipy.constants import physical_constants as cst
from warnings import warn
from math import pi
from . import elements
from ..lattice import checktype, AtWarning, AtError
from ..physics import linopt, ohmi_envelope

__all__ = ['Lattice']

TWO_PI_ERROR = 1.E-4


class Lattice(list):
    """Lattice object
    An AT lattice is a sequence of AT elements. This object accepts extended indexing"""

    _translate = dict(Energy='energy', Periodicity='periodicity', FamName='name')
    _ignore = set(('PassMethod', 'Length'))

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
            # Keep all attributes
            attrs = descr.__dict__
        else:
            params = [elem for elem in descr if isinstance(elem, elements.RingParam)]
            descr = [elem for elem in descr if not isinstance(elem, elements.RingParam)]
            radon = False
            for elem in descr:
                if elem.PassMethod.endswith('RadPass') or elem.PassMethod.endswith('CavityPass'):
                    radon = True
                    break
            # Initialize attributes to blank
            attrs = dict(name='', energy=None, periodicity=None, _radiation_on=radon)
            if len(params) > 0:
                # Set all RingParam attributes to the lattice object
                attrs.update(
                    (self._translate.get(key, key.lower()), value) for (key, value) in params[0].__dict__.items()
                    if key not in self._ignore)
        attrs.update(kwargs)

        super(Lattice, self).__init__(descr)

        if attrs['energy'] is None:
            # Guess energy from the Energy attribute of the elements
            # Look first in cavities
            energies = [elem.Energy for elem in descr if isinstance(elem, elements.RFCavity)]
            if len(energies) == 0:
                # Then look in all elements
                energies = [elem.Energy for elem in descr if hasattr(elem, 'Energy')]
            if len(energies) == 0:
                raise AtError('Energy not defined')
            energy = max(energies)
            if min(energies) < energy:
                warn(AtWarning('Inconsistent energy values, "Energy" set to {0}'.format(energy)))
            attrs['energy'] = energy

        if attrs['periodicity'] is None:
            # Guess periodicity from the bending angle of the superperiod
            theta = [elem.BendingAngle for elem in descr if isinstance(elem, elements.Dipole)]
            try:
                nbp = 2.0 * pi / sum(theta)
            except ZeroDivisionError:
                warn(AtWarning('No bending in the cell, set "Periodicity" to 1'))
                attrs['periodicity'] = 1
            else:
                periodicity = int(round(nbp))
                if abs(periodicity - nbp) > TWO_PI_ERROR:
                    warn(AtWarning('non-integer number of cells: {0} -> {1}'.format(nbp, periodicity)))
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
        """Energy loss per turn [eV]

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

    def linopt(self, *args, **kwargs):
        """
        Perform linear analysis of a lattice

        lindata0, tune, chrom[, lindata] = Lattice.linopt(dp[, refpts])

        PARAMETERS
            dp              momentum deviation. Defaults to 0
            refpts          Optional: elements at which data is returned. It can be
                            1) an integer (0 indicating the first element)
                            2) a list of integers
                            3) a numpy array of booleans as long as ring where
                               selected elements are true
                            Defaults to None

        KEYWORDS
            orbit           avoids looking for the colsed orbit if is already known ((6,) array)
            get_chrom=False compute dispersion and chromaticities. Needs computing the optics
                            at 2 different momentum deviations around the central one.
                            Defaults to False
            keep_lattice    Assume no lattice change since the previous tracking.
                            Defaults to False
            ddp=1.0E-8      momentum deviation used for computation of chromaticities and dispersion
            coupled=True    if False, simplify the calculations by assuming no H/V coupling

        OUTPUT
            lindata0        linear optics data at the entrance/end of the ring
            tune            [tune_A, tune_B], linear tunes for the two normal modes of linear motion [1]
            chrom           [ksi_A , ksi_B], vector of chromaticities ksi = d(nu)/(dP/P).
                            Only computed if 'get_chrom' is True
            lindata         Only returned if refpts is not None:
                            linear optics at the points refered to by refpts

            lindata is a structured array with fields:
            idx             element index in the ring                           (nrefs,)
            s_pos           longitudinal position [m]                           (nrefs,)
            closed_orbit    closed orbit vector with                            (nrefs, 6)
            dispersion      dispersion vector.                                  (nrefs, 4)
                            Only computed if 'get_chrom' is True                (nrefs, 4)
            m44             4x4 transfer matrix M from the beginning of ring    (nrefs, 4, 4)
                            to the entrance of the element [2]
            A               (2, 2) matrix A in [3]                              (nrefs, 2, 2)
            B               (2, 2) matrix B in [3]                              (nrefs, 2, 2)
            C               (2, 2) matrix C in [3]                              (nrefs, 2, 2)
            gamma           gamma parameter of the transformation to eigenmodes (nrefs,)
            mu              [mux, muy], A and B betatron phase (modulo 2*pi)    (nrefs, 2)
            beta            [betax, betay] vector                               (nrefs, 2)
            alpha           [alphax, alphay] vector                             (nrefs, 2)
            All values are given at the entrance of each element specified in refpts.

        REFERENCES
            [1] D.Edwars,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
            [2] E.Courant, H.Snyder
            [3] D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams, vol.2 (1999)

        See also get_twiss

        """
        if self._radiation_on:
            raise AtError('linopt needs no radiation in the lattice')
        return linopt(self, *args, **kwargs)

    def ohmi_envelope(self, *args, **kwargs):
        """
        Calculate the equilibrium beam envelope in a
        circular accelerator using Ohmi's beam envelope formalism [1]

        emit0, mode_emit, damping_rates, tunes[, emit] = ohmi_envelope([refpts])

        PARAMETERS
            refpts              elements at which data is returned. It can be
                                1) an integer (0 indicating the first element)
                                2) a list of integers
                                3) a numpy array of booleans as long as ring where selected elements are true
                                Defaults to None

        KEYWORDS
            orbit=None          Avoids looking for the colsed orbit if is already known ((6,) array)
            keep_lattice=False  Assume no lattice change since the previous tracking.
                                Defaults to False

        OUTPUT
            emit0               emittance data at the start/end of the ring
            beamdata            beam parameters at the start of the ring
            emit                Only returned if refpts is not None:
                                emittance data at the points refered to by refpts

            emit is a structured array with fields:
            R66                 (6, 6) equilibrium envelope matrix R
            R44                 (4, 4) betatron emittance matrix (dpp = 0)
            T66                 (6, 6) transfer matrix from the start of the ring
            orbit6              (6,) closed orbit
            emitXY              betatron emittance projected on xxp and yyp
            emitXYZ             6x6 emittance projected on xxp, yyp, ldp

            beamdata is a named tuple with attributes:
            tunes               tunes of the 3 normal modes
            damping_rates       damping rates of the 3 normal modes
            mode_matrices       R-matrices of the 3 normal modes
            mode_emittances     equilibrium emittances of the 3 normal modes

        REFERENCES
            [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
        """
        if not self._radiation_on:
            raise AtError('ohmi_envelope needs radiation in the lattice')
        return ohmi_envelope(self, *args, **kwargs)
