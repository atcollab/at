"""
Module to define common elements used in AT.

Each element has a default PassMethod attribute for which it should have the
appropriate attributes.  If a different PassMethod is set, it is the caller's
responsibility to ensure that the appropriate attributes are present.
"""
import re
import numpy
import copy
from inspect import getmembers, isdatadescriptor


def _array(value, shape=(-1,), dtype=numpy.float64):
    # Ensure proper ordering(F) and alignment(A) for "C" access in integrators
    return numpy.require(value, dtype=dtype, requirements=['F', 'A']).reshape(
        shape, order='F')


def _array66(value):
    return _array(value, shape=(6, 6))


def _nop(value):
    return value


class Element(object):
    """Base of pyat elements"""

    REQUIRED_ATTRIBUTES = ['FamName']
    _conversions = dict(FamName=str, PassMethod=str, Length=float,
                        R1=_array66, R2=_array66,
                        T1=lambda v: _array(v, (6,)),
                        T2=lambda v: _array(v, (6,)),
                        RApertures=lambda v: _array(v, (4,)),
                        EApertures=lambda v: _array(v, (2,)),
                        KickAngle=lambda v: _array(v, (2,)),
                        PolynomB=_array, PolynomA=_array,
                        BendingAngle=float,
                        MaxOrder=int, NumIntSteps=int,
                        Energy=float,
                        )

    _entrance_fields = ['T1', 'R1']
    _exit_fields = ['T2', 'R2']

    def __init__(self, family_name, **kwargs):
        self.FamName = family_name
        self.Length = kwargs.pop('Length', 0.0)
        self.PassMethod = kwargs.pop('PassMethod', 'IdentityPass')
        self.update(kwargs)

    def __setattr__(self, key, value):
        try:
            super(Element, self).__setattr__(
                key, self._conversions.get(key, _nop)(value))
        except Exception as exc:
            exc.args = ('In element {0}, parameter {1}: {2}'.format(
                self.FamName, key, exc),)
            raise

    def __str__(self):
        first3 = ['FamName', 'Length', 'PassMethod']
        attrs = dict(self.items())
        keywords = ['\t{0} : {1!s}'.format(k, attrs.pop(k)) for k in first3]
        keywords += ['\t{0} : {1!s}'.format(k, v) for k, v in attrs.items()]
        return '\n'.join((self.__class__.__name__ + ':', '\n'.join(keywords)))

    def __repr__(self):
        attrs = dict(self.items())
        arguments = [attrs.pop(k, getattr(self, k)) for k in
                     self.REQUIRED_ATTRIBUTES]
        defelem = self.__class__(*arguments)
        keywords = ['{0!r}'.format(arg) for arg in arguments]
        keywords += ['{0}={1!r}'.format(k, v) for k, v in sorted(attrs.items())
                     if not numpy.array_equal(v, getattr(defelem, k, None))]
        args = re.sub('\n\s*', ' ', ', '.join(keywords))
        return '{0}({1})'.format(self.__class__.__name__, args)

    def equals(self, other):
        """Whether an element is equivalent to another.

        This implementation was found to be too slow for the generic
        __eq__ method when comparing lattices.

        """
        return repr(self) == repr(other)

    def divide(self, frac):
        """split the element in len(frac) pieces whose length
        is frac[i]*self.Length

        arguments:
            frac            length of each slice expressed as a fraction of the
                            initial length. sum(frac) may differ from 1.

        Return a list of elements equivalent to the original.

        Example:

        >>> Drift('dr', 0.5).divide([0.2, 0.6, 0.2])
        [Drift('dr', 0.1), Drift('dr', 0.3), Drift('dr', 0.1)]
        """
        # Bx default, the element is indivisible
        return [self]

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
        """Return a shallow copy of the element"""
        return copy.copy(self)

    def deepcopy(self):
        """Return a deep copy of the element"""
        return copy.deepcopy(self)

    def items(self):
        """Iterates through the data members including slots and properties"""
        # Get attributes
        for k, v in vars(self).items():
            yield k, v
        # Get slots and properties
        for k, v in getmembers(self.__class__, isdatadescriptor):
            if not k.startswith('_'):
                yield k, getattr(self, k)


class LongElement(Element):
    """pyAT long element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length']

    def __init__(self, family_name, length, *args, **kwargs):
        kwargs.setdefault('Length', length)
        # Ancestor may be either Element of ThinMultipole
        # noinspection PyArgumentList
        super(LongElement, self).__init__(family_name, *args, **kwargs)

    def _part(self, fr, sumfr):
        pp = self.copy()
        pp.Length = fr * self.Length
        if hasattr(self, 'KickAngle'):
            pp.KickAngle = fr / sumfr * self.KickAngle
        return pp

    def divide(self, frac):
        """split the element in len(frac) pieces whose length
        is frac[i]*self.Length

        arguments:
            frac            length of each slice expressed as a fraction of the
                            initial length. sum(frac) may differ from 1.

        Return a list of elements equivalent to the original.

        Example:

        >>> Drift('dr', 0.5).divide([0.2, 0.6, 0.2])
        [Drift('dr', 0.1), Drift('dr', 0.3), Drift('dr', 0.1)]
        """

        def popattr(element, attr):
            val = getattr(element, attr)
            delattr(element, attr)
            return attr, val

        frac = numpy.asarray(frac, dtype=float)
        el = self.copy()
        # Remove entrance and exit attributes
        fin = dict(popattr(el, key) for key in vars(self) if
                   key in self._entrance_fields)
        fout = dict(popattr(el, key) for key in vars(self) if
                    key in self._exit_fields)
        # Split element
        element_list = [el._part(f, numpy.sum(frac)) for f in frac]
        # Restore entrance and exit attributes
        for key, value in fin.items():
            setattr(element_list[0], key, value)
        for key, value in fout.items():
            setattr(element_list[-1], key, value)
        return element_list


class Marker(Element):
    """pyAT marker element"""


class Monitor(Element):
    """pyAT monitor element"""


class Aperture(Element):
    """pyAT aperture element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Limits']
    _conversions = dict(Element._conversions, Limits=lambda v: _array(v, (4,)))

    def __init__(self, family_name, limits, **kwargs):
        kwargs.setdefault('PassMethod', 'AperturePass')
        super(Aperture, self).__init__(family_name, Limits=limits, **kwargs)


class Drift(LongElement):
    """pyAT drift space element"""

    def __init__(self, family_name, length, **kwargs):
        """Drift(FamName, Length, **keywords)
        """
        kwargs.setdefault('PassMethod', 'DriftPass')
        super(Drift, self).__init__(family_name, length, **kwargs)

    def insert(self, insert_list):
        """insert elements inside a drift

        arguments:
            insert_list: iterable, each item of insert_list is itself an
                         iterable with 2 objects:
                             1. the location where the center of the element
                                will be inserted, given as a fraction of the
                                Drift length.
                             2. an element to be inserted at that location. If
                                None, the drift will be divided but no element
                                will be inserted.

        Return a list of elements.

        Drifts with negative lengths may be generated if necessary.

        Examples:

        >>> Drift('dr', 2.0).insert(((0.25, None), (0.75, None)))
        [Drift('dr', 0.5), Drift('dr', 1.0), Drift('dr', 0.5)]

        >>> Drift('dr', 2.0).insert(((0.0, Marker('m1')), (0.5, Marker('m2'))))
        [Marker('m1'), Drift('dr', 1.0), Marker('m2'), Drift('dr', 1.0)]

        >>> Drift('dr', 2.0).insert(((0.5, Quadrupole('qp', 0.4, 0.0)),))
        [Drift('dr', 0.8), Quadrupole('qp', 0.4), Drift('dr', 0.8)]
        """
        frac, elements = zip(*insert_list)
        lg = [0.0 if el is None else el.Length for el in elements]
        fr = numpy.asarray(frac, dtype=float)
        lg = 0.5 * numpy.asarray(lg, dtype=float) / self.Length
        drfrac = numpy.hstack((fr - lg, 1.0)) - numpy.hstack((0.0, fr + lg))
        long_elems = (drfrac != 0.0)
        drifts = numpy.ndarray((len(drfrac),), dtype='O')
        drifts[long_elems] = self.divide(drfrac[long_elems])
        line = [None] * (len(drifts) + len(elements))
        line[::2] = drifts
        line[1::2] = elements
        return [el for el in line if el is not None]


class ThinMultipole(Element):
    """pyAT thin multipole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['PolynomA',
                                                         'PolynomB']

    def __init__(self, family_name, poly_a, poly_b, **kwargs):
        """ThinMultipole(FamName, PolynomA, PolynomB, **keywords)

        Available keywords:
        MaxOrder        Number of desired multipoles. Default: highest index of
                        non-zero polynomial coefficients
        """
        def getpol(poly):
            nonzero = numpy.flatnonzero(poly != 0.0)
            return poly, len(poly), nonzero[-1] if len(nonzero) > 0 else -1

        def lengthen(poly, dl):
            if dl > 0:
                return numpy.concatenate((poly, numpy.zeros(dl)))
            else:
                return poly

        # Remove MaxOrder, PolynomA and PolynomB
        poly_a, len_a, ord_a = getpol(_array(kwargs.pop('PolynomA', poly_a)))
        poly_b, len_b, ord_b = getpol(_array(kwargs.pop('PolynomB', poly_b)))
        deforder = max(getattr(self, 'DefaultOrder', 0), ord_a, ord_b)
        maxorder = kwargs.pop('MaxOrder', deforder)
        len_ab = max(maxorder + 1, len_a, len_b)
        kwargs.setdefault('PassMethod', 'ThinMPolePass')
        super(ThinMultipole, self).__init__(family_name, **kwargs)
        # Set MaxOrder while PolynomA and PolynomB are not set yet
        super(ThinMultipole, self).__setattr__('MaxOrder', maxorder)
        # Adjust polynom lengths and set them
        self.PolynomA = lengthen(poly_a, len_ab - len_a)
        self.PolynomB = lengthen(poly_b, len_ab - len_b)

    def __setattr__(self, key, value):
        """Check the compatibility of MaxOrder, PolynomA and PolynomB"""
        polys = ('PolynomA', 'PolynomB')
        if key in polys:
            value = _array(value)
            lmin = getattr(self, 'MaxOrder')
            if not len(value) > lmin:
                raise ValueError(
                    'Length of {0} must be larger than {1}'.format(key, lmin))
        elif key == 'MaxOrder':
            value = int(value)
            lmax = min(len(getattr(self, k)) for k in polys)
            if not value < lmax:
                raise ValueError(
                    'MaxOrder must be smaller than {0}'.format(lmax))

        super(ThinMultipole, self).__setattr__(key, value)


class Multipole(LongElement, ThinMultipole):
    """pyAT multipole element"""
    REQUIRED_ATTRIBUTES = LongElement.REQUIRED_ATTRIBUTES + ['PolynomA',
                                                             'PolynomB']
    _conversions = dict(ThinMultipole._conversions, K=float, H=float)

    def __init__(self, family_name, length, poly_a, poly_b, **kwargs):
        """Multipole(FamName, Length, PolynomA, PolynomB, **keywords)

        Available keywords:
        MaxOrder        Number of desired multipoles. Default: highest index of
                        non-zero polynomial coefficients
        NumIntSteps     Number of integration steps (default: 10)
        KickAngle       Correction deviation angles (H, V)
        """
        kwargs.setdefault('PassMethod', 'StrMPoleSymplectic4Pass')
        kwargs.setdefault('NumIntSteps', 10)
        super(Multipole, self).__init__(family_name, length,
                                        poly_a, poly_b, **kwargs)


class Dipole(Multipole):
    """pyAT dipole element"""
    REQUIRED_ATTRIBUTES = LongElement.REQUIRED_ATTRIBUTES + ['BendingAngle',
                                                             'K']
    _conversions = dict(Multipole._conversions, EntranceAngle=float,
                        ExitAngle=float,
                        FringeInt1=float, FringeInt2=float,
                        FringeQuadEntrance=int, FringeQuadExit=int,
                        FringeBendEntrance=int, FringeBendExit=int)

    _entrance_fields = Multipole._entrance_fields + ['EntranceAngle',
                                                     'FringeInt1',
                                                     'FringeBendEntrance'
                                                     'FringeQuadEntrance']
    _exit_fields = Multipole._exit_fields + ['ExitAngle',
                                             'FringeInt2',
                                             'FringeBendExit',
                                             'FringeQuadExit']

    DefaultOrder = 0

    def __init__(self, family_name, length, bending_angle=0.0, k=0.0,
                 **kwargs):
        """Dipole(FamName, Length, bending_angle, Strength=0, **keywords)

        Available keywords:
        EntranceAngle   entrance angle (default 0.0)
        ExitAngle       exit angle (default 0.0)
        PolynomB        straight multipoles
        PolynomA        skew multipoles
        MaxOrder        Number of desired multipoles
        NumIntSteps     Number of integration steps (default: 10)
        FullGap         Magnet full gap
        FringeInt1      Fringe field extension
        FringeInt2
        FringeBendEntrance  1: legacy version Brown First Order (default)
                            2: SOLEIL close to second order of Brown
                            3: THOMX
        FringeBendExit
        FringeQuadEntrance  0: no fringe fiels effect (default)
                            1: Lee-Whiting's thin lens limit formula
                            2: elegant-like
        FringeQuadExit
        fringeIntM0     Integrals for FringeQuad method 2
        fringeIntP0
        KickAngle       Correction deviation angles (H, V)
        """
        poly_b = kwargs.pop('PolynomB', numpy.array([0, k]))
        kwargs.setdefault('BendingAngle', bending_angle)
        kwargs.setdefault('EntranceAngle', 0.0)
        kwargs.setdefault('ExitAngle', 0.0)
        kwargs.setdefault('PassMethod', 'BendLinearPass')
        super(Dipole, self).__init__(family_name, length, [], poly_b, **kwargs)

    def _part(self, fr, sumfr):
        pp = super(Dipole, self)._part(fr, sumfr)
        pp.BendingAngle = fr / sumfr * self.BendingAngle
        pp.EntranceAngle = 0.0
        pp.ExitAngle = 0.0
        return pp

    # noinspection PyPep8Naming
    @property
    def K(self):
        return self.PolynomB[1]

    # noinspection PyPep8Naming
    @K.setter
    def K(self, strength):
        self.PolynomB[1] = strength


# Bend is a synonym of Dipole.
Bend = Dipole


class Quadrupole(Multipole):
    """pyAT quadrupole element"""
    REQUIRED_ATTRIBUTES = LongElement.REQUIRED_ATTRIBUTES + ['K']
    _conversions = dict(Multipole._conversions, FringeQuadEntrance=int,
                        FringeQuadExit=int)

    _entrance_fields = Multipole._entrance_fields + ['FringeQuadEntrance']
    _exit_fields = Multipole._exit_fields + ['FringeQuadExit']

    DefaultOrder = 1

    def __init__(self, family_name, length, k=0.0, **kwargs):
        """Quadrupole(FamName, Length, Strength=0, **keywords)

        Available keywords:
        PolynomB        straight multipoles
        PolynomA        skew multipoles
        MaxOrder        Number of desired multipoles
        NumIntSteps     Number of integration steps (default: 10)
        FringeQuadEntrance  0: no fringe fiels effect (default)
                            1: Lee-Whiting's thin lens limit formula
                            2: elegant-like
        FringeQuadExit
        fringeIntM0     Integrals for FringeQuad method 2
        fringeIntP0
        KickAngle       Correction deviation angles (H, V)
        """
        poly_b = kwargs.pop('PolynomB', numpy.array([0, k]))
        kwargs.setdefault('PassMethod', 'QuadLinearPass')
        super(Quadrupole, self).__init__(family_name, length, [], poly_b,
                                         **kwargs)

    # noinspection PyPep8Naming
    @property
    def K(self):
        return self.PolynomB[1]

    # noinspection PyPep8Naming
    @K.setter
    def K(self, strength):
        self.PolynomB[1] = strength


class Sextupole(Multipole):
    """pyAT sextupole element"""
    REQUIRED_ATTRIBUTES = LongElement.REQUIRED_ATTRIBUTES + ['H']

    DefaultOrder = 2

    def __init__(self, family_name, length, h=0.0, **kwargs):
        """Sextupole(FamName, Length, Strength=0, **keywords)

        Available keywords:
        PolynomB        straight multipoles
        PolynomA        skew multipoles
        MaxOrder        Number of desired multipoles
        NumIntSteps     Number of integration steps (default: 10)
        KickAngle       Correction deviation angles (H, V)
        """
        poly_b = kwargs.pop('PolynomB', [0, 0, h])
        kwargs.setdefault('PassMethod', 'StrMPoleSymplectic4Pass')
        super(Sextupole, self).__init__(family_name, length, [], poly_b,
                                        **kwargs)

    # noinspection PyPep8Naming
    @property
    def H(self):
        return self.PolynomB[2]

    # noinspection PyPep8Naming
    @H.setter
    def H(self, strength):
        self.PolynomB[2] = strength


class Octupole(Multipole):
    """pyAT octupole element, with no changes from multipole at present"""
    REQUIRED_ATTRIBUTES = Multipole.REQUIRED_ATTRIBUTES

    DefaultOrder = 3


class RFCavity(LongElement):
    """pyAT RF cavity element"""
    REQUIRED_ATTRIBUTES = LongElement.REQUIRED_ATTRIBUTES + ['Voltage',
                                                             'Frequency',
                                                             'HarmNumber',
                                                             'Energy']
    _conversions = dict(LongElement._conversions,
                        Voltage=float, Frequency=float,
                        HarmNumber=int, TimeLag=float)

    def __init__(self, family_name, length, voltage, frequency,
                 harmonic_number, energy, **kwargs):
        """
        Available keywords:
        TimeLag   time lag with respect to the reference particle
        """
        kwargs.setdefault('TimeLag', 0.0)
        kwargs.setdefault('PassMethod', 'CavityPass')
        super(RFCavity, self).__init__(family_name, length,
                                       Voltage=voltage,
                                       Frequency=frequency,
                                       HarmNumber=harmonic_number,
                                       Energy=energy, **kwargs)

    def _part(self, fr, sumfr):
        pp = super(RFCavity, self)._part(fr, sumfr)
        pp.Voltage = fr * self.Voltage
        return pp


class M66(Element):
    """Linear (6, 6) transfer matrix"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES
    _conversions = dict(Element._conversions, M66=_array66)

    def __init__(self, family_name, m66=None, **kwargs):
        if m66 is None:
            m66 = numpy.identity(6)
        kwargs.setdefault('PassMethod', 'Matrix66Pass')
        super(M66, self).__init__(family_name, M66=m66, **kwargs)


class Corrector(LongElement):
    """pyAT corrector element"""
    REQUIRED_ATTRIBUTES = LongElement.REQUIRED_ATTRIBUTES + ['KickAngle']

    def __init__(self, family_name, length, kick_angle, **kwargs):
        kwargs.setdefault('PassMethod', 'CorrectorPass')
        super(Corrector, self).__init__(family_name, length,
                                        KickAngle=kick_angle, **kwargs)


class Wiggler(LongElement):
    """pyAT wiggler element

    See atwiggler.m
    """
    REQUIRED_ATTRIBUTES = LongElement.REQUIRED_ATTRIBUTES + ['Lw', 'Bmax',
                                                             'Energy']
    _conversions = dict(Element._conversions, Lw=float, Bmax=float,
                        Energy=float,
                        Bx=lambda v: _array(v, (6, -1)),
                        By=lambda v: _array(v, (6, -1)),
                        Nstep=int, Nmeth=int, NHharm=int, NVharm=int)

    def __init__(self, family_name, length, wiggle_period, b_max, energy,
                 Nstep=5, Nmeth=4, By=(1, 1, 0, 1, 1, 0), Bx=(), **kwargs):
        """
        Args:
            length: total length of the wiggler
            wiggle_period: length must be a multiple of this
            b_max: peak wiggler field [Tesla]
            energy: beam energy [eV]

        Available keywords:
            Nstep: number of integration steps.
            Nmeth: symplectic integration order: 2 or 4
            Bx: harmonics for horizontal wiggler: (6,nHharm) array-like object
            By: harmonics for vertical wiggler (6,nHharm) array-like object

        """
        kwargs.setdefault('PassMethod', 'GWigSymplecticPass')
        n_wiggles = length / wiggle_period
        if abs(round(n_wiggles) - n_wiggles) > 1e-6:
            raise ValueError("Wiggler: length / wiggle_period is not an "
                             "integer. ({0}/{1}={2})".format(length,
                                                             wiggle_period,
                                                             n_wiggles))
        super(Wiggler, self).__init__(family_name, length, Lw=wiggle_period,
                                      Bmax=b_max, Nstep=Nstep, Nmeth=Nmeth,
                                      By=By, Bx=Bx, Energy=energy, **kwargs)

        for i, b in enumerate(self.By.T):
            dk = abs(b[3]**2 - b[4]**2 - b[2]**2) / abs(b[4])
            if dk > 1e-6:
                raise ValueError("Wiggler(H): kx^2 + kz^2 -ky^2 !=0, i = "
                                 "{0}".format(i))

        for i, b in enumerate(self.Bx.T):
            dk = abs(b[2]**2 - b[4]**2 - b[3]**2) / abs(b[4])
            if dk > 1e-6:
                raise ValueError("Wiggler(V): ky^2 + kz^2 -kx^2 !=0, i = "
                                 "{0}".format(i))

        self.NHharm = self.By.shape[1]
        self.NVharm = self.Bx.shape[1]


CLASS_MAP = dict((k, v) for k, v in locals().items()
                 if isinstance(v, type) and issubclass(v, Element))
