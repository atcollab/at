"""
Module to define common elements used in AT.

Each element has a default PassMethod attribute for which it should have the
appropriate attributes.  If a different PassMethod is set, it is the caller's
responsibility to ensure that the appropriate attributes are present.
"""
import numpy
import itertools


class Element(object):
    REQUIRED_ATTRIBUTES = ['FamName']
    FLOAT_ARRAYS = ['R1', 'R2', 'T1', 'T2']

    def __init__(self, family_name, length=0.0, **kwargs):
        self.FamName = family_name
        self.Length = float(length)
        self.PassMethod = kwargs.pop('PassMethod', 'IdentityPass')
        for field in Element.FLOAT_ARRAYS:
            if field in kwargs:
                kwargs[field] = numpy.ascontiguousarray(kwargs[field],
                                                        dtype=numpy.float64)
        for k in kwargs:
            setattr(self, k, kwargs[k])

    def __str__(self):
        keywords = ('{0:>16} : {1!r}'.format(k, v) for k, v in self.__dict__.items())
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
        return '{0}({1})'.format(self.__class__.__name__, ','.join(itertools.chain(arguments, keywords)))


class Marker(Element):
    """pyAT marker element"""

    def __init__(self, family_name, **kwargs):
        super(Marker, self).__init__(family_name, kwargs.pop('Length', 0.0), **kwargs)


class Monitor(Element):
    """pyAT monitor element"""

    def __init__(self, family_name, **kwargs):
        super(Monitor, self).__init__(family_name, kwargs.pop('Length', 0.0), **kwargs)


class Aperture(Element):
    """pyAT aperture element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Limits']

    def __init__(self, family_name, limits, **kwargs):
        assert len(limits) == 4
        kwargs['Limits'] = numpy.array(limits, dtype=numpy.float64)
        kwargs.setdefault('PassMethod', 'AperturePass')
        super(Aperture, self).__init__(family_name, kwargs.pop('Length', 0.0), **kwargs)


class Drift(Element):
    """pyAT drift space element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length']

    def __init__(self, family_name, length, **kwargs):
        """Drift(FamName, Length, **keywords)
        """
        kwargs.setdefault('PassMethod', 'DriftPass')
        super(Drift, self).__init__(family_name, length, **kwargs)


class ThinMultipole(Element):
    """pyAT thin multipole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['PolynomA',
                                                         'PolynomB']

    def __init__(self, family_name, poly_a, poly_b, **kwargs):
        """ThinMultipole(FamName, PolynomA, PolynomB, **keywords)

        Available keywords:
        'MaxOrder'      Number of desired multipoles
        """
        poly_a = numpy.array(kwargs.pop('PolynomA', poly_a), dtype=numpy.float64)
        poly_b = numpy.array(kwargs.pop('PolynomB', poly_b), dtype=numpy.float64)
        poly_size = max(kwargs.get('MaxOrder', 0) + 1, len(poly_a), len(poly_b))
        kwargs['PolynomA'] = numpy.concatenate((poly_a, numpy.zeros(poly_size - len(poly_a))))
        kwargs['PolynomB'] = numpy.concatenate((poly_b, numpy.zeros(poly_size - len(poly_b))))
        kwargs['MaxOrder'] = int(kwargs.pop('MaxOrder', poly_size - 1))
        kwargs.setdefault('PassMethod', 'ThinMPolePass')
        super(ThinMultipole, self).__init__(family_name, kwargs.pop('Length', 0.0), **kwargs)


class Multipole(ThinMultipole):
    """pyAT multipole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length',
                                                         'PolynomA',
                                                         'PolynomB']

    def __init__(self, family_name, length, poly_a, poly_b, **kwargs):
        """Multipole(FamName, Length, PolynomA, PolynomB, **keywords)

        Available keywords:
        'MaxOrder'      Number of desired multipoles
        'NumIntSteps'   Number of integration steps (default: 10)
        """
        kwargs['NumIntSteps'] = int(kwargs.pop('NumIntSteps', 10))
        kwargs['Length'] = length
        kwargs.setdefault('PassMethod', 'StrMPoleSymplectic4Pass')
        super(Multipole, self).__init__(family_name, poly_a, poly_b, **kwargs)


class Dipole(Multipole):
    """pyAT dipole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length',
                                                         'BendingAngle']

    def __init__(self, family_name, length, bending_angle, k=0.0, **kwargs):
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
        kwargs['BendingAngle'] = float(bending_angle)
        kwargs['EntranceAngle'] = float(kwargs.pop('EntranceAngle', 0.0))
        kwargs['ExitAngle'] = float(kwargs.pop('ExitAngle', 0.0))
        kwargs.setdefault('PassMethod', 'BendLinearPass')
        super(Dipole, self).__init__(family_name, length, [], poly_b, **kwargs)


# Bend is a synonym of Dipole.
Bend = Dipole


class Quadrupole(Multipole):
    """pyAT quadrupole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length']

    def __init__(self, family_name, length, k=0.0, **kwargs):
        """Quadrupole(FamName, Length, Strength=0, **keywords)

        Available keywords:
        'PolynomB'      straight multipoles
        'PolynomA'      skew multipoles
        'MaxOrder'      Number of desired multipoles
        'NumIntSteps'   Number of integration steps (default: 10)
        """
        poly_b = kwargs.pop('PolynomB', numpy.array([0, k]))
        kwargs.setdefault('PassMethod', 'QuadLinearPass')
        super(Quadrupole, self).__init__(family_name, length, [], poly_b, **kwargs)


class Sextupole(Multipole):
    """pyAT sextupole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length']

    def __init__(self, family_name, length, h=0.0, **kwargs):
        """Sextupole(FamName, Length, Strength=0, **keywords)

        Available keywords:
        'PolynomB'  straight multipoles
        'PolynomA'  skew multipoles
        'MaxOrder'  Number of desired multipoles
        'NumIntSteps'   Number of integration steps (default: 10)
        """
        poly_b = kwargs.pop('PolynomB', [0, 0, h])
        kwargs.setdefault('PassMethod', 'StrMPoleSymplectic4Pass')
        super(Sextupole, self).__init__(family_name, length, [], poly_b, **kwargs)


class Octupole(Multipole):
    """pyAT octupole element, with no changes from multipole at present"""
    pass


class RFCavity(Element):
    """pyAT RF cavity element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length',
                                                         'Voltage',
                                                         'Frequency',
                                                         'HarmNumber',
                                                         'Energy']

    def __init__(self, family_name, length, voltage, frequency, harmonic_number, energy, **kwargs):
        """
        Available keywords:
        'TimeLag'   time lag with respect to the reference particle
        """
        kwargs.setdefault('Voltage', float(voltage))
        kwargs.setdefault('Frequency', float(frequency))
        kwargs.setdefault('HarmNumber', int(harmonic_number))
        kwargs.setdefault('Energy', float(energy))
        kwargs.setdefault('TimeLag', 0.0)
        kwargs.setdefault('PassMethod', 'CavityPass')
        super(RFCavity, self).__init__(family_name, length, **kwargs)


class RingParam(Element):
    """pyAT RingParam element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Energy',
                                                         'Periodicity']

    def __init__(self, family_name, energy, nb_periods, **kwargs):
        kwargs.setdefault('Energy', float(energy))
        kwargs.setdefault('Periodicity', int(nb_periods))
        kwargs.setdefault('PassMethod', 'IdentityPass')
        super(RingParam, self).__init__(family_name, kwargs.pop('Length', 0.0), **kwargs)
