"""
Module to define common elements used in AT.

Each element has a default PassMethod attribute for which it should have the
appropriate attributes.  If a different PassMethod is set, it is the caller's
responsibility to ensure that the appropriate attributes are present.
"""
import numpy
import itertools


def _array(value, shape=(-1,), dtype=numpy.float64):
    # Ensure proper ordering(F) and alignment(A) for "C" access in integrators
    return numpy.require(value, dtype=dtype, requirements=['F', 'A']).reshape(shape)


def _array66(value):
    return _array(value, shape=(6, 6))


def _float(value):
    return float(value)


def _int(value):
    return int(value)


def _nop(value):
    return value


class Element(object):
    REQUIRED_ATTRIBUTES = ['FamName']
    CONVERSIONS = dict(R1=_array66, R2=_array66, T1=lambda v: _array(v, (6,)), T2=lambda v: _array(v, (6,)),
                       RApertures=lambda v: _array(v, (4,)), EApertures=lambda v: _array(v, (2,)))

    def __init__(self, family_name, Length=0.0, PassMethod='IdentityPass', **kwargs):
        self.FamName = family_name
        self.Length = float(Length)
        self.PassMethod = PassMethod

        for (key, value) in kwargs.items():
            try:
                setattr(self, key, self.CONVERSIONS.get(key, _nop)(value))
            except Exception as exc:
                exc.args = ('In element {0}, parameter {1}: {2}'.format(family_name, key, exc),)
                raise

    def __str__(self):
        first3 = ['FamName', 'Length', 'PassMethod']
        keywords = ['{0} : {1!s}'.format(k, self.__dict__[k]) for k in first3]
        keywords = keywords + ['{0} : {1!s}'.format(k, v) for k, v in
                               self.__dict__.items() if k not in first3]
        return '\n'.join((self.__class__.__name__ + ':', '\n'.join(keywords)))

    def __repr__(self):
        def differ(v1, v2):
            if isinstance(v1, numpy.ndarray):
                return not numpy.array_equal(v1, v2)
            else:
                return v1 != v2

        defelem = self.__class__(*(getattr(self, k) for k in self.REQUIRED_ATTRIBUTES))
        arguments = ('{0!r}'.format(getattr(self, k)) for k in self.REQUIRED_ATTRIBUTES)
        keywords = ('{0}={1!r}'.format(k, v) for k, v in self.__dict__.items() if differ(v, getattr(defelem, k, None)))
        return '{0}({1})'.format(self.__class__.__name__, ', '.join(itertools.chain(arguments, keywords)))


class Marker(Element):
    """pyAT marker element"""

    def __init__(self, family_name, **kwargs):
        super(Marker, self).__init__(family_name, **kwargs)


class Monitor(Element):
    """pyAT monitor element"""

    def __init__(self, family_name, **kwargs):
        super(Monitor, self).__init__(family_name, **kwargs)


class Aperture(Element):
    """pyAT aperture element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Limits']
    CONVERSIONS = dict(Element.CONVERSIONS, Limits=lambda v: _array(v, (4,)))

    def __init__(self, family_name, limits, **kwargs):
        kwargs.setdefault('PassMethod', 'AperturePass')
        super(Aperture, self).__init__(family_name, Limits=limits, **kwargs)


class Drift(Element):
    """pyAT drift space element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length']

    def __init__(self, family_name, length, **kwargs):
        """Drift(FamName, Length, **keywords)
        """
        kwargs.setdefault('PassMethod', 'DriftPass')
        super(Drift, self).__init__(family_name, Length=kwargs.pop('Length', length), **kwargs)


class ThinMultipole(Element):
    """pyAT thin multipole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['PolynomA',
                                                         'PolynomB']
    CONVERSIONS = dict(Element.CONVERSIONS, BendingAngle=_float, MaxOrder=_int,
                       PolynomB=_array, PolynomA=_array)

    def __init__(self, family_name, poly_a, poly_b, MaxOrder=0, **kwargs):
        """ThinMultipole(FamName, PolynomA, PolynomB, **keywords)

        Available keywords:
        'MaxOrder'      Number of desired multipoles
        """
        poly_a = numpy.array(kwargs.pop('PolynomA', poly_a), dtype=numpy.float64)
        poly_b = numpy.array(kwargs.pop('PolynomB', poly_b), dtype=numpy.float64)
        poly_size = max(MaxOrder + 1, len(poly_a), len(poly_b))
        poly_a = numpy.concatenate((poly_a, numpy.zeros(poly_size - len(poly_a))))
        poly_b = numpy.concatenate((poly_b, numpy.zeros(poly_size - len(poly_b))))
        kwargs.setdefault('PassMethod', 'ThinMPolePass')
        super(ThinMultipole, self).__init__(family_name, PolynomB=poly_b, PolynomA=poly_a,
                                            MaxOrder=MaxOrder, **kwargs)


class Multipole(ThinMultipole):
    """pyAT multipole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length',
                                                         'PolynomA',
                                                         'PolynomB']
    CONVERSIONS = dict(ThinMultipole.CONVERSIONS, NumIntSteps=_int, KickAngle=_array)

    def __init__(self, family_name, length, poly_a, poly_b, NumIntSteps=10, **kwargs):
        """Multipole(FamName, Length, PolynomA, PolynomB, **keywords)

        Available keywords:
        'MaxOrder'      Number of desired multipoles
        'NumIntSteps'   Number of integration steps (default: 10)
        """
        kwargs.setdefault('PassMethod', 'StrMPoleSymplectic4Pass')
        super(Multipole, self).__init__(family_name, poly_a, poly_b, Length=kwargs.pop('Length', length),
                                        NumIntSteps=NumIntSteps, **kwargs)


class Dipole(Multipole):
    """pyAT dipole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length', 'BendingAngle']
    CONVERSIONS = dict(Multipole.CONVERSIONS, EntranceAngle=_float, ExitAngle=_float,
                       FringeQuadEntrance=_int, FringeQuadExit=_int,
                       FringeBendEntrance=_int, FringeBendExit=_int)

    def __init__(self, family_name, length, BendingAngle, k=0.0, EntranceAngle=0.0, ExitAngle=0.0, **kwargs):
        """Dipole(FamName, Length, BendingAngle, Strength=0, **keywords)

        Available keywords:
        'EntranceAngle' entrance angle (default 0.0)
        'ExitAngle'     exit angle (default 0.0)
        'PolynomB'      straight multipoles
        'PolynomA'      skew multipoles
        'MaxOrder'      Number of desired multipoles
        'NumIntSteps'   Number of integration steps (default: 10)
        """
        poly_b = kwargs.pop('PolynomB', numpy.array([0, k]))
        kwargs.setdefault('PassMethod', 'BendLinearPass')
        super(Dipole, self).__init__(family_name, length, [], poly_b, BendingAngle=BendingAngle,
                                     EntranceAngle=EntranceAngle, ExitAngle=ExitAngle, **kwargs)


# Bend is a synonym of Dipole.
Bend = Dipole


class Quadrupole(Multipole):
    """pyAT quadrupole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length']
    CONVERSIONS = dict(Multipole.CONVERSIONS, FringeQuadEntrance=_int, FringeQuadExit=_int)

    def __init__(self, family_name, length, k=0.0, MaxOrder=1, **kwargs):
        """Quadrupole(FamName, Length, Strength=0, **keywords)

        Available keywords:
        'PolynomB'      straight multipoles
        'PolynomA'      skew multipoles
        'MaxOrder'      Number of desired multipoles
        'NumIntSteps'   Number of integration steps (default: 10)
        """
        poly_b = kwargs.pop('PolynomB', numpy.array([0, k]))
        kwargs.setdefault('PassMethod', 'QuadLinearPass')
        super(Quadrupole, self).__init__(family_name, length, [], poly_b, MaxOrder=MaxOrder, **kwargs)


class Sextupole(Multipole):
    """pyAT sextupole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length']

    def __init__(self, family_name, length, h=0.0, MaxOrder=2, **kwargs):
        """Sextupole(FamName, Length, Strength=0, **keywords)

        Available keywords:
        'PolynomB'  straight multipoles
        'PolynomA'  skew multipoles
        'MaxOrder'  Number of desired multipoles
        'NumIntSteps'   Number of integration steps (default: 10)
        """
        poly_b = kwargs.pop('PolynomB', [0, 0, h])
        kwargs.setdefault('PassMethod', 'StrMPoleSymplectic4Pass')
        super(Sextupole, self).__init__(family_name, length, [], poly_b, MaxOrder=MaxOrder, **kwargs)


class Octupole(Multipole):
    """pyAT octupole element, with no changes from multipole at present"""
    REQUIRED_ATTRIBUTES = Multipole.REQUIRED_ATTRIBUTES


class RFCavity(Element):
    """pyAT RF cavity element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length', 'Voltage', 'Frequency', 'HarmNumber', 'Energy']
    CONVERSIONS = dict(Element.CONVERSIONS, Voltage=_float, Frequency=_float, HarmNumber=_int, Energy=_float,
                       TimeLag=_float)

    def __init__(self, family_name, length, voltage, frequency, harmonic_number, energy, TimeLag=0.0, **kwargs):
        """
        Available keywords:
        'TimeLag'   time lag with respect to the reference particle
        """
        kwargs.setdefault('PassMethod', 'CavityPass')
        super(RFCavity, self).__init__(family_name, length, Voltage=voltage, Frequency=frequency,
                                       HarmNumber=harmonic_number, Energy=energy, TimeLag=TimeLag, **kwargs)


class RingParam(Element):
    """pyAT RingParam element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Energy']
    CONVERSIONS = dict(Element.CONVERSIONS, Energy=_float, Periodicity=_int)

    def __init__(self, family_name, energy, Periodicity=1, **kwargs):
        kwargs.setdefault('PassMethod', 'IdentityPass')
        super(RingParam, self).__init__(family_name, Energy=energy, Periodicity=Periodicity, **kwargs)


class M66(Element):
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES
    CONVERSIONS = dict(Element.CONVERSIONS, M66=_array66)

    def __init__(self, family_name, m66=None, **kwargs):
        if m66 is None:
            m66 = numpy.identity(6)
        kwargs.setdefault('PassMethod', 'Matrix66Pass')
        super(M66, self).__init__(family_name, M66=m66, **kwargs)


class Corrector(Element):
    """pyAT corrector element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length', 'KickAngle']
    CONVERSIONS = dict(Element.CONVERSIONS, KickAngle=_array)

    def __init__(self, family_name, length, kick_angle, **kwargs):
        kwargs.setdefault('PassMethod', 'CorrectorPass')
        super(Corrector, self).__init__(family_name, kwargs.pop('Length', length),
                                        KickAngle=kick_angle, **kwargs)
