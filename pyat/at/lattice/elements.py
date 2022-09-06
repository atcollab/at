"""
Module to define common elements used in AT.

Each element has a default PassMethod attribute for which it should have the
appropriate attributes.  If a different PassMethod is set, it is the caller's
responsibility to ensure that the appropriate attributes are present.
"""
import abc
import re
import numpy
from copy import copy, deepcopy
from abc import ABC
from typing import Optional, Generator, Tuple, List, Iterable


def _array(value, shape=(-1,), dtype=numpy.float64):
    # Ensure proper ordering(F) and alignment(A) for "C" access in integrators
    return numpy.require(value, dtype=dtype, requirements=['F', 'A']).reshape(
        shape, order='F')


def _array66(value):
    return _array(value, shape=(6, 6))


def _nop(value):
    return value


class LongtMotion(ABC):
    """Abstract Base class for all Element classes whose instances may modify
    the particle momentum

    Allows to identify elements potentially inducing longitudinal motion.

    Subclasses of :py:class:`LongtMotion` must provide two methods for
    enabling longitudinal motion:

    * ``_get_longt_motion(self)`` must return the activation state,
    * ``set_longt_motion(self, enable, new_pass=None, copy=False, **kwargs)``
      must enable or disable longitudinal motion.
    """
    @abc.abstractmethod
    def _get_longt_motion(self):
        return False

    # noinspection PyShadowingNames
    @abc.abstractmethod
    def set_longt_motion(self, enable, new_pass=None, copy=False, **kwargs):
        """Enable/Disable longitudinal motion

        Parameters:
            enable:     ``True`` for enabling, ``False`` for disabling
            new_pass:   New PassMethod:

              * ``None``: makes no change,
              * ``'auto'``: Uses the default conversion,
              * Anything else is used as the new PassMethod.
            copy:       If True, returns a modified copy of the element,
              otherwise modifies the element in-place
        """
        # noinspection PyUnresolvedReferences
        if new_pass is None or new_pass == self.PassMethod:
            return self if copy else None
        if copy:
            newelem = deepcopy(self)
            newelem.PassMethod = new_pass
            return newelem
        # noinspection PyAttributeOutsideInit
        self.PassMethod = new_pass


# noinspection PyUnresolvedReferences
class _DictLongtMotion(LongtMotion):
    # noinspection PyShadowingNames
    """Mixin class for elements implementing a 'default_pass' class attribute

    :py:class:`DictLongtMotion` provides:

    * a :py:meth:`set_longt_motion` method setting the PassMethod according
      to the ``default_pass`` dictionary.
    * a :py:obj:`.longt_motion` property set to ``True`` when the PassMethod
      is ``default_pass[True]``

    The class must have a ``default_pass`` class attribute, a dictionary
    such that:

    * ``default_pass[False]`` is the PassMethod when radiation is turned
      OFF,
    * ``default_pass[True]`` is the default PassMethod when radiation is
      turned ON.

    The :py:class:`DictLongtMotion` class must be set as the first base class.

    Example:

        >>> class QuantumDiffusion(_DictLongtMotion, Element):
        ...
        ...     default_pass = {False: 'IdentityPass', True: 'QuantDiffPass'}

        Defines a class such that :py:meth:`set_longt_motion` will select
        ``'IdentityPass'`` or ``'IdentityPass'``.
        """
    def _get_longt_motion(self):
        return self.PassMethod != self.default_pass[False]

    # noinspection PyShadowingNames
    def set_longt_motion(self, enable, new_pass=None, **kwargs):
        if new_pass == 'auto':
            new_pass = self.default_pass[enable]
        return super().set_longt_motion(enable, new_pass=new_pass, **kwargs)


# noinspection PyUnresolvedReferences
class _Radiative(LongtMotion):
    # noinspection PyShadowingNames
    r"""Mixin class for radiating elements

    :py:class:`_Radiative` implements the mechanism for converting the pass
    methods of radiating elements. It provides:

    * a :py:meth:`set_longt_motion` method setting the PassMethod
      according to the following rule:

      * ``enable == True``: replace "\*Pass" by "\*RadPass"
      * ``enable == False``: replace "\*RadPass" by "\*Pass"
    * a :py:obj:`.longt_motion` property set to true when the PassMethod
      ends with "RadPass"

    The :py:class:`_Radiative` class must be set as the first base class.

    Example:
        >>> class Multipole(_Radiative, LongElement, ThinMultipole):

        Defines a class where :py:meth:`set_longt_motion` will convert the
        PassMethod according to the \*Pass or \*RadPass suffix.
    """

    def _get_longt_motion(self):
        return self.PassMethod.endswith('RadPass')

    def _autopass(self, enable):
        rad = self.longt_motion
        if enable and not rad:
            return ''.join((self.PassMethod[:-4], 'RadPass'))
        elif not enable and rad:
            return ''.join((self.PassMethod[:-7], 'Pass'))
        else:
            return None

    # noinspection PyTypeChecker,PyShadowingNames
    def set_longt_motion(self, enable, new_pass=None, copy=False, **kwargs):
        if new_pass == 'auto':
            new_pass = self._autopass(enable)
        if new_pass is None or new_pass == self.PassMethod:
            return self if copy else None
        if enable:
            def setpass(el):
                el.PassMethod = new_pass
                el.Energy = kwargs['energy']
        else:
            def setpass(el):
                el.PassMethod = new_pass
                try:
                    del el.Energy
                except AttributeError:
                    pass
        if copy:
            newelem = deepcopy(self)
            setpass(newelem)
            return newelem
        setpass(self)


class Radiative(_Radiative):
    # noinspection PyShadowingNames
    r"""Mixin class for default radiating elements (:py:class:`.Dipole`,
    :py:class:`.Quadrupole`, :py:class:`.Wiggler`)

    :py:class:`Radiative` is a base class for the subset of radiative elements
    considered as the ones to be turned on by default: :py:class:`.Dipole`,
    :py:class:`.Quadrupole` and :py:class:`.Wiggler`, excluding the higher
    order multipoles.

    :py:class:`Radiative` inherits from :py:class:`_Radiative` and does not
    add any new functionality. Its purpose is to identify the default set of
    radiating elements.

    Example:
        >>> class Dipole(Radiative, Multipole):

        Defines a class belonging to the default radiating elements. It
        converts the PassMethod according to the "\*Pass" or "\*RadPass"
        suffix.
    """
    pass


class Collective(_DictLongtMotion):
    """Mixin class for elements representing collective effects

    Derived classes will automatically set the :py:obj:`is_collective`
    property when the element is active.

    The class must have a ``default_pass`` class attribute, a dictionary such
    that:

    * ``default_pass[False]`` is the PassMethod when collective effects
      are turned OFF,
    * ``default_pass[True]`` is the default PassMethod when collective effects
      are turned ON.

    The :py:class:`Collective` class must be set as the first base class.

    Example:
        >>> class WakeElement(Collective, Element):
        ...
        ...     default_pass = {False: 'IdentityPass', True: 'WakeFieldPass'}

        Defines a class where the :py:obj:`is_collective` property is handled
    """

    def _get_collective(self):
        # noinspection PyUnresolvedReferences
        return self.PassMethod != self.default_pass[False]


class Element(object):
    """Base class for AT elements"""

    _BUILD_ATTRIBUTES = ['FamName']
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

    def __init__(self, family_name: str, **kwargs):
        """
        Parameters:
            family_name:    Name of the element

        All keywords will be set as attributes of the element
        """

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
        return '\n'.join((type(self).__name__ + ':', '\n'.join(keywords)))

    def __repr__(self):
        attrs = dict(self.items())
        arguments = [attrs.pop(k, getattr(self, k)) for k in
                     self._BUILD_ATTRIBUTES]
        defelem = self.__class__(*arguments)
        keywords = ['{0!r}'.format(arg) for arg in arguments]
        keywords += ['{0}={1!r}'.format(k, v) for k, v in sorted(attrs.items())
                     if not numpy.array_equal(v, getattr(defelem, k, None))]
        args = re.sub(r'\n\s*', ' ', ', '.join(keywords))
        return '{0}({1})'.format(self.__class__.__name__, args)

    def equals(self, other) -> bool:
        """Whether an element is equivalent to another.

        This implementation was found to be too slow for the generic
        __eq__ method when comparing lattices.
        """
        return repr(self) == repr(other)

    def divide(self, frac) -> List["Element"]:
        """split the element in len(frac) pieces whose length
        is frac[i]*self.Length

        Parameters:
            frac:           length of each slice expressed as a fraction of the
              initial length. ``sum(frac)`` may differ from 1.

        Returns:
            elem_list:  a list of elements equivalent to the original.

        Example:

            >>> Drift('dr', 0.5).divide([0.2, 0.6, 0.2])
            [Drift('dr', 0.1), Drift('dr', 0.3), Drift('dr', 0.1)]
        """
        # Bx default, the element is indivisible
        return [self]

    def update(self, *args, **kwargs):
        """
        update(**kwargs)
        update(mapping, **kwargs)
        update(iterable, **kwargs)
        Update the element attributes with the given arguments
        """
        attrs = dict(*args, **kwargs)
        for (key, value) in attrs.items():
            setattr(self, key, value)

    def copy(self) -> "Element":
        """Return a shallow copy of the element"""
        return copy(self)

    def deepcopy(self) -> "Element":
        """Return a deep copy of the element"""
        return deepcopy(self)

    def items(self) -> Generator[Tuple, None, None]:
        """Iterates through the data members"""
        for k, v in vars(self).items():
            yield k, v

    def is_compatible(self, other) -> bool:
        """Checks if another Element can be merged"""
        return False

    def merge(self, other) -> None:
        """Merge another element"""
        if not self.is_compatible(other):
            badname = getattr(other, 'FamName', type(other))
            raise TypeError('Cannot merge {0} and {1}'.format(self.FamName,
                                                              badname))

    # noinspection PyMethodMayBeStatic
    def _get_longt_motion(self):
        return False

    # noinspection PyMethodMayBeStatic
    def _get_collective(self):
        return False

    @property
    def longt_motion(self):
        """``True`` if longitudinal motion is affected by the element"""
        return self._get_longt_motion()

    @property
    def is_collective(self):
        """``True`` if the element involves collective effects"""
        return self._get_collective()


class LongElement(Element):
    """Base class for long elements
    """
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ['Length']

    def __init__(self, family_name: str, length: float, *args, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]

        Other arguments and keywords are given to the base class
        """
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

    def divide(self, frac) -> List[Element]:
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

    def is_compatible(self, other) -> bool:
        return type(other) is type(self) and \
               self.PassMethod == other.PassMethod

    def merge(self, other) -> None:
        super().merge(other)
        self.Length += other.Length


class Marker(Element):
    """Marker element"""


class Monitor(Element):
    """Monitor element"""


class Aperture(Element):
    """Aperture element"""
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ['Limits']
    _conversions = dict(Element._conversions, Limits=lambda v: _array(v, (4,)))

    def __init__(self, family_name, limits, **kwargs):
        """
        Args:
            family_name:    Name of the element
            limits:         (4,) array of physical aperture:
              [xmin, xmax, zmin, zmax] [m]

        Default PassMethod: ``AperturePass``
        """
        kwargs.setdefault('PassMethod', 'AperturePass')
        super(Aperture, self).__init__(family_name, Limits=limits, **kwargs)


class Drift(LongElement):
    """Drift space element"""

    def __init__(self, family_name: str, length: float, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]

        Default PassMethod: ``DriftPass``
        """
        kwargs.setdefault('PassMethod', 'DriftPass')
        super(Drift, self).__init__(family_name, length, **kwargs)

    def insert(self,
               insert_list: Iterable[Tuple[float, Optional[Element]]]) \
            -> List[Element]:
        """insert elements inside a drift

        Arguments:
            insert_list: iterable, each item of insert_list is itself an
              iterable with 2 objects:

              1. the location where the center of the element
                 will be inserted, given as a fraction of the Drift length.
              2. an element to be inserted at that location. If ``None``,
                 the drift will be divided but no element will be inserted.

        Returns:
             elem_list: a list of elements.

        Drifts with negative lengths may be generated if necessary.

        Examples:

            >>> Drift('dr', 2.0).insert(((0.25, None), (0.75, None)))
            [Drift('dr', 0.5), Drift('dr', 1.0), Drift('dr', 0.5)]

            >>> Drift('dr', 2.0).insert(((0.0, Marker('m1')),
            ... (0.5, Marker('m2'))))
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
        nline = len(drifts) + len(elements)
        line = [None] * nline  # type: List[Optional[Element]]
        line[::2] = drifts
        line[1::2] = elements
        return [el for el in line if el is not None]


class Collimator(Drift):
    """Collimator element"""
    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ['RApertures']

    def __init__(self, family_name: str, length: float, limits, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            limits:         (4,) array of physical aperture:
              [xmin, xmax, zmin, zmax] [m]

        Default PassMethod: ``DriftPass``
        """
        super(Collimator, self).__init__(family_name, length,
                                         RApertures=limits, **kwargs)


class ThinMultipole(Element):
    """Thin multipole element"""
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ['PolynomA',
                                                     'PolynomB']

    def __init__(self, family_name: str, poly_a, poly_b, **kwargs):
        """
        Args:
            family_name:    Name of the element
            poly_a:         Array of normal multipole components
            poly_b:         Array of skew multipole components

        Keyword arguments:
            MaxOrder:   Number of desired multipoles. Default: highest index of
                        non-zero polynomial coefficients

        Default PassMethod: ``ThinMPolePass``
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
        kwargs.setdefault('PassMethod', 'ThinMPolePass')
        super(ThinMultipole, self).__init__(family_name, **kwargs)
        # Set MaxOrder while PolynomA and PolynomB are not set yet
        super(ThinMultipole, self).__setattr__('MaxOrder', maxorder)
        # Adjust polynom lengths and set them
        len_ab = max(self.MaxOrder + 1, len_a, len_b)
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


class Multipole(_Radiative, LongElement, ThinMultipole):
    """Multipole element"""
    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ['PolynomA',
                                                         'PolynomB']
    _conversions = dict(ThinMultipole._conversions, K=float, H=float)

    def __init__(self, family_name: str, length: float, poly_a, poly_b,
                 **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            poly_a:         Array of normal multipole components
            poly_b:         Array of skew multipole components

        Keyword arguments:
            MaxOrder:   Number of desired multipoles. Default: highest index of
              non-zero polynomial coefficients
            NumIntSteps: Number of integration steps (default: 10)
            KickAngle:  Correction deviation angles (H, V)

        Default PassMethod: ``StrMPoleSymplectic4Pass``
        """
        kwargs.setdefault('PassMethod', 'StrMPoleSymplectic4Pass')
        kwargs.setdefault('NumIntSteps', 10)
        super(Multipole, self).__init__(family_name, length,
                                        poly_a, poly_b, **kwargs)

    def is_compatible(self, other) -> bool:
        if super().is_compatible(other) and \
                self.MaxOrder == other.MaxOrder:
            for i in range(self.MaxOrder + 1):
                if self.PolynomB[i] != other.PolynomB[i]:
                    return False
                if self.PolynomA[i] != other.PolynomA[i]:
                    return False
            return True
        else:
            return False

    # noinspection PyPep8Naming
    @property
    def K(self) -> float:
        """Focusing strength [mˆ-2]"""
        return 0.0 if len(self.PolynomB) < 2 else self.PolynomB[1]

    # noinspection PyPep8Naming
    @K.setter
    def K(self, strength: float):
        self.PolynomB[1] = strength

    # noinspection PyPep8Naming
    @property
    def H(self) -> float:
        """Sextupolar strength"""
        return 0.0 if len(self.PolynomB) < 3 else self.PolynomB[2]

    # noinspection PyPep8Naming
    @H.setter
    def H(self, strength):
        self.PolynomB[2] = strength


class Dipole(Radiative, Multipole):
    """Dipole element"""
    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ['BendingAngle',
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

    def __init__(self, family_name: str, length: float,
                 bending_angle: Optional[float] = 0.0, k: float = 0.0,
                 **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            bending_angle:  Bending angle [rd]
            poly_a:         Array of normal multipole components
            poly_b:         Array of skew multipole components
            k=0:            Field index

        Keyword arguments:
            EntranceAngle=0.0:  entrance angle
            ExitAngle=0.0:      exit angle
            PolynomB:           straight multipoles
            PolynomA:           skew multipoles
            MaxOrder:           Number of desired multipoles
            NumIntSt=10:        Number of integration steps
            FullGap:            Magnet full gap
            FringeInt1:         Fringe field extension
            FringeInt2:
            FringeBendEntrance: 1: legacy version Brown First Order (default)

              2: SOLEIL close to second order of Brown

              3: THOMX
            FringeBendExit:     See ``FringeBendEntrance``
            FringeQuadEntrance: 0: no fringe field effect (default)

              1: Lee-Whiting's thin lens limit formula

              2: elegant-like
            FringeQuadExit:     See ``FringeQuadEntrance``
            fringeIntM0:        Integrals for FringeQuad method 2
            fringeIntP0:
            KickAngle:          Correction deviation angles (H, V)

        Default PassMethod: ``BndMPoleSymplectic4Pass``
        """
        poly_b = kwargs.pop('PolynomB', numpy.array([0, k]))
        kwargs.setdefault('BendingAngle', bending_angle)
        kwargs.setdefault('EntranceAngle', 0.0)
        kwargs.setdefault('ExitAngle', 0.0)
        kwargs.setdefault('PassMethod', 'BndMPoleSymplectic4Pass')
        super(Dipole, self).__init__(family_name, length, [], poly_b, **kwargs)

    def items(self) -> Generator[Tuple, None, None]:
        yield from super().items()
        yield 'K', self.K

    def _part(self, fr, sumfr):
        pp = super(Dipole, self)._part(fr, sumfr)
        pp.BendingAngle = fr / sumfr * self.BendingAngle
        pp.EntranceAngle = 0.0
        pp.ExitAngle = 0.0
        return pp

    def is_compatible(self, other) -> bool:
        def invrho(dip: Dipole):
            return dip.BendingAngle / dip.Length

        return (super().is_compatible(other) and
                self.ExitAngle == -other.EntranceAngle and
                abs(invrho(self) - invrho(other)) <= 1.e-6)

    def merge(self, other) -> None:
        super().merge(other)
        # noinspection PyAttributeOutsideInit
        self.ExitAngle = other.ExitAngle
        self.BendingAngle += other.BendingAngle


# Bend is a synonym of Dipole.
Bend = Dipole


class Quadrupole(Radiative, Multipole):
    """Quadrupole element"""
    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ['K']
    _conversions = dict(Multipole._conversions, FringeQuadEntrance=int,
                        FringeQuadExit=int)

    _entrance_fields = Multipole._entrance_fields + ['FringeQuadEntrance']
    _exit_fields = Multipole._exit_fields + ['FringeQuadExit']

    DefaultOrder = 1

    def __init__(self, family_name: str, length: float,
                 k: Optional[float] = 0.0, **kwargs):
        """Quadrupole(FamName, Length, Strength=0, **keywords)

        Args:
            family_name:    Name of the element
            length:         Element length [m]
            k:              strength [mˆ-2]

        Keyword Arguments:
            PolynomB:           straight multipoles
            PolynomA:           skew multipoles
            MaxOrder:           Number of desired multipoles
            NumIntSteps=10:     Number of integration steps
            FringeQuadEntrance: 0: no fringe field effect (default)

              1: Lee-Whiting's thin lens limit formula

              2: elegant-like
            FringeQuadExit:     See ``FringeQuadEntrance``
            fringeIntM0:        Integrals for FringeQuad method 2
            fringeIntP0:
            KickAngle:          Correction deviation angles (H, V)

        Default PassMethod: ``StrMPoleSymplectic4Pass``
        """
        poly_b = kwargs.pop('PolynomB', numpy.array([0, k]))
        kwargs.setdefault('PassMethod', 'StrMPoleSymplectic4Pass')
        super(Quadrupole, self).__init__(family_name, length, [], poly_b,
                                         **kwargs)

    def items(self) -> Generator[Tuple, None, None]:
        yield from super().items()
        yield 'K', self.K


class Sextupole(Multipole):
    """Sextupole element"""
    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ['H']

    DefaultOrder = 2

    def __init__(self, family_name: str, length: float,
                 h: Optional[float] = 0.0, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            h:              strength [mˆ-3]

        Keyword Arguments:
            PolynomB:           straight multipoles
            PolynomA:           skew multipoles
            MaxOrder:           Number of desired multipoles
            NumIntSteps=10:     Number of integration steps
            KickAngle:          Correction deviation angles (H, V)

        Default PassMethod: ``StrMPoleSymplectic4Pass``
        """
        poly_b = kwargs.pop('PolynomB', [0, 0, h])
        kwargs.setdefault('PassMethod', 'StrMPoleSymplectic4Pass')
        super(Sextupole, self).__init__(family_name, length, [], poly_b,
                                        **kwargs)

    def items(self) -> Generator[Tuple, None, None]:
        yield from super().items()
        yield 'H', self.H


class Octupole(Multipole):
    """Octupole element, with no changes from multipole at present"""
    _BUILD_ATTRIBUTES = Multipole._BUILD_ATTRIBUTES

    DefaultOrder = 3


class RFCavity(LongtMotion, LongElement):
    """RF cavity element"""
    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ['Voltage',
                                                         'Frequency',
                                                         'HarmNumber',
                                                         'Energy']
    default_pass = {False: 'DriftPass', True: 'RFCavityPass'}
    _conversions = dict(LongElement._conversions,
                        Voltage=float, Frequency=float,
                        HarmNumber=int, TimeLag=float)

    def __init__(self, family_name: str, length: float, voltage: float,
                 frequency: float, harmonic_number: int, energy: float,
                 **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            voltage:        RF voltage [V]
            frequency:      RF frequency [Hz]
            harmonic_number:
            energy:         ring energy [eV]

        Keyword Arguments:
            TimeLag=0:      Cavity time lag

        Default PassMethod: ``RFCavityPass``
        """
        kwargs.setdefault('TimeLag', 0.0)
        kwargs.setdefault('PassMethod', self.default_pass[True])
        super(RFCavity, self).__init__(family_name, length,
                                       Voltage=voltage,
                                       Frequency=frequency,
                                       HarmNumber=harmonic_number,
                                       Energy=energy, **kwargs)

    def _part(self, fr, sumfr):
        pp = super(RFCavity, self)._part(fr, sumfr)
        pp.Voltage = fr * self.Voltage
        return pp

    def is_compatible(self, other) -> bool:
        return (super().is_compatible(other) and
                self.Frequency == other.Frequency and
                self.TimeLag == other.TimeLag)

    def merge(self, other) -> None:
        super().merge(other)
        self.Voltage += other.Voltage

    def _get_longt_motion(self):
        return self.PassMethod.endswith('CavityPass')

    # noinspection PyShadowingNames
    def set_longt_motion(self, enable, new_pass=None, **kwargs):
        if new_pass == 'auto':
            new_pass = (self.default_pass[True] if enable else
                        ('IdentityPass' if self.Length == 0 else 'DriftPass'))
        return super().set_longt_motion(enable, new_pass=new_pass, **kwargs)


class M66(Element):
    """Linear (6, 6) transfer matrix"""
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES
    _conversions = dict(Element._conversions, M66=_array66)

    def __init__(self, family_name: str, m66=None, **kwargs):
        """
        Args:
            family_name:    Name of the element
            m66:            Transfer matrix. Default: Identity matrix

        Default PassMethod: ``Matrix66Pass``
       """
        if m66 is None:
            m66 = numpy.identity(6)
        kwargs.setdefault('PassMethod', 'Matrix66Pass')
        super(M66, self).__init__(family_name, M66=m66, **kwargs)


class Corrector(LongElement):
    """Corrector element"""
    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ['KickAngle']

    def __init__(self, family_name: str, length: float, kick_angle, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            KickAngle:      Correction deviation angles (H, V)

        Default PassMethod: ``CorrectorPass``
        """
        kwargs.setdefault('PassMethod', 'CorrectorPass')
        super(Corrector, self).__init__(family_name, length,
                                        KickAngle=kick_angle, **kwargs)


class Wiggler(Radiative, LongElement):
    """Wiggler element

    See atwiggler.m
    """
    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ['Lw', 'Bmax',
                                                         'Energy']
    _conversions = dict(Element._conversions, Lw=float, Bmax=float,
                        Energy=float,
                        Bx=lambda v: _array(v, (6, -1)),
                        By=lambda v: _array(v, (6, -1)),
                        Nstep=int, Nmeth=int, NHharm=int, NVharm=int)

    # noinspection PyPep8Naming
    def __init__(self, family_name: str, length: float, wiggle_period: float,
                 b_max: float, energy: float, Nstep: Optional[int] = 5,
                 Nmeth: Optional[int] = 4,
                 By=(1, 1, 0, 1, 1, 0), Bx=(), **kwargs):
        """
        Args:
            length:         total length of the wiggler
            wiggle_period:  length must be a multiple of this
            b_max:          peak wiggler field [Tesla]
            energy:         beam energy [eV]
            Nstep:          number of integration steps.
            Nmeth:          symplectic integration order: 2 or 4
            Bx:             harmonics for horizontal wiggler: (6, nHharm)
                              array-like object
            By:             harmonics for vertical wiggler (6, nHharm)
                              array-like object

        Default PassMethod: ``GWigSymplecticPass``
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
            dk = abs(b[3] ** 2 - b[4] ** 2 - b[2] ** 2) / abs(b[4])
            if dk > 1e-6:
                raise ValueError("Wiggler(H): kx^2 + kz^2 -ky^2 !=0, i = "
                                 "{0}".format(i))

        for i, b in enumerate(self.Bx.T):
            dk = abs(b[2] ** 2 - b[4] ** 2 - b[3] ** 2) / abs(b[4])
            if dk > 1e-6:
                raise ValueError("Wiggler(V): ky^2 + kz^2 -kx^2 !=0, i = "
                                 "{0}".format(i))

        self.NHharm = self.By.shape[1]
        self.NVharm = self.Bx.shape[1]


class QuantumDiffusion(_DictLongtMotion, Element):
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ['Lmatp']
    default_pass = {False: 'IdentityPass', True: 'QuantDiffPass'}
    _conversions = dict(Element._conversions, Lmatp=_array66)

    def __init__(self, family_name: str, lmatp: numpy.ndarray, **kwargs):
        """Quantum diffusion element

        Args:
            family_name:    Name of the element
            lmatp      :    Diffusion matrix for generation (see
                               ``at.physics.radiation.gen_quandiff_elem``)

        Default PassMethod: ``QuantDiffPass``
        """
        kwargs.setdefault('PassMethod', self.default_pass[True])
        super().__init__(family_name, Lmatp=lmatp, **kwargs)


def build_class_map():  # Missing class aliases (Bend)
    global CLASS_MAP

    def subclasses_recursive(cl):
        direct = cl.__subclasses__()
        indirect = []
        for subclass in direct:
            indirect.extend(subclasses_recursive(subclass))
        return frozenset([cl] + direct + indirect)

    cls_list = subclasses_recursive(Element)
    CLASS_MAP = {cls.__name__: cls for cls in cls_list}


def get_class_map():
    return CLASS_MAP


# build_class_map()

CLASS_MAP = dict((k, v) for k, v in locals().items()
                 if isinstance(v, type) and issubclass(v, Element))
