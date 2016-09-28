import numpy


class Element(object):
    REQUIRED_ATTRIBUTES = ['FamName']

    def __init__(self, family_name, length=0.0, **kwargs):
        self.FamName = family_name
        self.Length = length
        self.PassMethod = kwargs.pop('PassMethod', 'IdentityPass')
        for k in kwargs:
            setattr(self, k, kwargs[k])


class Marker(Element):
    """pyAT marker element"""

    def __init__(self, family_name, **kwargs):
        super(Marker, self).__init__(family_name, 0.0, **kwargs)


class Monitor(Element):
    """pyAT monitor element"""

    def __init__(self, family_name, **kwargs):
        super(Monitor, self).__init__(family_name, 0.0, **kwargs)


class Aperture(Element):
    """pyAT aperture element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Limits']

    def __init__(self, family_name, limits, **kwargs):
        assert len(limits) == 4
        kwargs['Limits'] = numpy.array(limits, dtype=numpy.float64)
        kwargs.setdefault('PassMethod', 'AperturePass')
        super(Aperture, self).__init__(family_name, 0.0, **kwargs)


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
        sz = max(kwargs.get('MaxOrder', 0) + 1, len(poly_a), len(poly_b))
        kwargs['PolynomA'] = numpy.concatenate((poly_a, numpy.zeros(sz - len(poly_a))))
        kwargs['PolynomB'] = numpy.concatenate((poly_b, numpy.zeros(sz - len(poly_b))))
        kwargs['MaxOrder'] = int(kwargs.pop('MaxOrder', sz - 1))
        lg = kwargs.pop('Length', 0.0)
        super(ThinMultipole, self).__init__(family_name, lg, **kwargs)


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
        super(Multipole, self).__init__(family_name, poly_a, poly_b, **kwargs)


class Dipole(Multipole):
    """pyAT dipole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length',
                                                         'BendingAngle']

    def __init__(self, family_name, length, bending_angle, K=0.0, **kwargs):
        """Dipole(FamName, Length, BendingAngle, Strength=0, **keywords)

        Available keywords:
        'EntranceAngle' entrance angle (default 0.0)
        'ExitAngle'     exit angle (default 0.0)
        'PolynomB'      straight multipoles
        'PolynomA'      skew multipoles
        'MaxOrder'      Number of desired multipoles
        'NumIntSteps'   Number of integration steps (default: 10)
        """
        poly_b = kwargs.pop('PolynomB', [0, K])
        kwargs.setdefault('EntranceAngle', 0.0)
        kwargs.setdefault('ExitAngle', 0.0)
        kwargs.setdefault('PassMethod', 'BndMPoleSymplectic4E2Pass')
        super(Dipole, self).__init__(family_name, length, [], poly_b, BendingAngle=bending_angle, **kwargs)


class Bend(Dipole):
    pass


class Quadrupole(Multipole):
    """pyAT quadrupole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length']

    def __init__(self, family_name, length, K=0.0, **kwargs):
        """Quadrupole(FamName, Length, Strength=0, **keywords)

        Available keywords:
        'PolynomB'      straight multipoles
        'PolynomA'      skew multipoles
        'MaxOrder'      Number of desired multipoles
        'NumIntSteps'   Number of integration steps (default: 10)
        """
        poly_b = kwargs.pop('PolynomB', [0, K])
        kwargs.setdefault('PassMethod', 'QuadLinearPass')
        super(Quadrupole, self).__init__(family_name, length, [], poly_b, **kwargs)


class Sextupole(Multipole):
    """pyAT sextupole element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Length']

    def __init__(self, family_name, length, H=0.0, **kwargs):
        """Sextupole(FamName, Length, Strength=0, **keywords)

        Available keywords:
        'PolynomB'  straight multipoles
        'PolynomA'  skew multipoles
        'MaxOrder'  Number of desired multipoles
        'NumIntSteps'   Number of integration steps (default: 10)
        """
        poly_b = kwargs.pop('PolynomB', [0, 0, H])
        kwargs.setdefault('PassMethod', 'StrMPoleSymplectic4Pass')
        super(Sextupole, self).__init__(family_name, length, [], poly_b, **kwargs)


class RFCavity(Element):
    """pyAT RF cavity element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Voltage',
                                                         'Frequency',
                                                         'HarmNumber',
                                                         'Energy',
                                                         'TimeLag']

    def __init__(self, family_name, length, voltage, frequency, harmonic_number, energy, **kwargs):
        """
        Available keywords:
        'TimeLag'   time lag with respect to the reference particle
        """
        kwargs.setdefault('Voltage', voltage)
        kwargs.setdefault('Frequency', frequency)
        kwargs.setdefault('HarmNumber', harmonic_number)
        kwargs.setdefault('Energy', energy)
        kwargs.setdefault('TimeLag', 0.0)
        kwargs.setdefault('PassMethod', 'CavityPass')
        super(RFCavity, self).__init__(family_name, length, **kwargs)


class RingParam(Element):
    """pyAT RingParam element"""
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES + ['Energy',
                                                         'Periodicity']

    def __init__(self, family_name, energy, nb_periods, **kwargs):
        kwargs.setdefault('Energy', energy)
        kwargs.setdefault('Periodicity', nb_periods)
        kwargs.setdefault('PassMethod', 'IdentityPass')
        super(RingParam, self).__init__(family_name, 0.0, **kwargs)
