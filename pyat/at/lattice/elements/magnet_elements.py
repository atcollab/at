""":py:class:`.Element` classes for magnets"""

from __future__ import annotations

__all__ = [
    "ThinMultipole",
    "Multipole",
    "Dipole",
    "Bend",
    "Quadrupole",
    "Sextupole",
    "Octupole",
    "Corrector",
    "Wiggler",
]

import warnings
from collections.abc import Generator, Callable
from typing import Any
import contextlib
from math import factorial

import numpy as np

from ..exceptions import AtError, AtWarning
from .conversions import _float, _array
from .abstract_elements import Radiative, _Radiative
from .element_object import Element
from .basic_elements import LongElement

# AtWarning from this module should always be issued (not only on the first occurrence)
warnings.filterwarnings("always", category=AtWarning, module=__name__)


def _k_property(attrname: str, order: int, doc: str):
    """Return a property for a field component in MAD convention."""

    def kget(mult: ThinMultipole) -> float:
        if order > mult.MaxOrder:
            return 0.0
        else:
            return factorial(order) * float(getattr(mult, attrname)[order])

    def kset(mult: ThinMultipole, value: float) -> None:
        getattr(mult, attrname)[order] = np.array(value) / factorial(order)

    return property(kget, kset, doc=doc)


def _s_property(attrname: str, doc: str):
    """Return a property for the strength of the main field component."""

    def kget(mult: ThinMultipole) -> float:
        order = getattr(mult, "DefaultOrder", None)
        if order is None:
            msg = "Undefined order\nYou need to set a 'DefaultOrder' attribute"
            raise AttributeError(msg)
        elif order > mult.MaxOrder:
            return 0.0
        else:
            return float(getattr(mult, attrname)[order])

    def kset(mult: ThinMultipole, value: float) -> None:
        order = getattr(mult, "DefaultOrder", None)
        if order is None:
            msg = "Undefined order\nYou need to set a 'DefaultOrder' attribute"
            raise AttributeError(msg)
        getattr(mult, attrname)[order] = np.array(value) / factorial(order)

    return property(kget, kset, doc=doc)


class _LinkedArray(np.ndarray):
    """Linked array.

    A numpy array triggering an action on any modification of its elements.
    """

    setitem: Callable[[int, Any], None]

    def __setitem__(self, key, value):
        # normal action
        super().__setitem__(key, value)
        with contextlib.suppress(AttributeError):
            # external action, if one is defined
            self.setitem(key, value)

    def __repr__(self) -> str:
        # Simulate a standard ndarray
        return repr(self.view(np.ndarray))


class ThinMultipole(Element):
    """Thin multipole element"""

    # Class attributes
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["PolynomA", "PolynomB"]
    _conversions = dict(Element._conversions, K=float, H=float)
    _stacklevel = 4  # Stacklevel for warnings

    # Instance attributes
    PolynomA: np.ndarray  #: Integrated skew strengths
    PolynomB: np.ndarray  #: Integrated straight strengths

    def __init__(self, family_name: str, poly_a, poly_b, **kwargs):
        """
        Args:
            family_name:    Name of the element
            poly_a:         Array of integrated skew multipole components
            poly_b:         Array of integrated normal multipole components

        Keyword arguments:
            MaxOrder:       Number of desired multipoles. Default: highest
              index of non-zero polynomial coefficients
            FieldScaling:   Scaling factor applied to the magnetic field
              (*PolynomA* and *PolynomB*)

        Default PassMethod: ``ThinMPolePass``
        """

        def getpol(poly):
            nonzero = np.flatnonzero(poly != 0.0)
            return len(poly), nonzero[-1] if len(nonzero) > 0 else -1

        def lengthen(poly, dl):
            if dl > 0:
                return np.concatenate((poly, np.zeros(dl)))
            else:
                return poly

        def make_same_length(arr1, arr2):
            diff = len(arr1) - len(arr2)
            if diff > 0:
                arr2 = lengthen(arr2, diff)
            elif diff < 0:
                arr1 = lengthen(arr1, -diff)
            return arr1, arr2

        def seterr(name, kname, kval, aname, aval):
            mess = (
                f"Element {name}: Conflicting element data, {kname!r} ({kval}) "
                f"in kwargs does not match positional argument {aname!r} ({aval})."
            )
            return AtError(mess)

        def setwarn(name, aname, kname):
            mess = (
                f"Element {name}: Duplicate element data, both positional "
                f"argument {aname!r} and {kname!r} in kwargs passed."
            )
            warnings.warn(AtWarning(mess), stacklevel=self._stacklevel)

        def check_polynom(keyname, arg):
            argvalue = self._conversions[keyname](arg)
            argname = f"poly_{keyname[-1].lower()}"
            if keyname in kwargs:
                kvalue = self._conversions[keyname](kwargs.pop(keyname))
                kvalue, argvalue = make_same_length(kvalue, argvalue)
                if issubclass(self.__class__, (Dipole, Quadrupole)):
                    if (
                        keyname == "PolynomB"
                        and argvalue[1] != 0.0
                        and (kvalue.size < 2 or argvalue[1] != kvalue[1])
                    ):
                        raise seterr(family_name, "PolynomB", kvalue, "k", argvalue[1])
                elif issubclass(self.__class__, Sextupole):
                    if (
                        keyname == "PolynomB"
                        and argvalue[2] != 0.0
                        and (kvalue.size < 3 or argvalue[2] != kvalue[2])
                    ):
                        raise seterr(family_name, "PolynomB", kvalue, "h", argvalue[2])
                elif np.any(argvalue) and not np.array_equiv(kvalue, argvalue):
                    raise seterr(family_name, keyname, kvalue, argname, argvalue)
                else:
                    setwarn(family_name, argname, keyname)
                return kvalue
            else:
                return argvalue

        def check_strength(keyname, index):
            if keyname in kwargs:
                k = self._conversions[keyname](kwargs.pop(keyname))
                vname = f"poly_b[{index}]"
                if len(prmpolb) > index:
                    if k != 0.0 and k != prmpolb[index]:
                        raise seterr(family_name, "K", k, vname, prmpolb[index])
                    else:
                        setwarn(family_name, "poly_b", "K")
                else:
                    raise seterr(family_name, "K", k, vname, 0.0)

        # To avoid potentially unintended behaviour, due to the constructor being
        # passed multiple definitions of the same thing, we do the following:
        # - If it is present in kwargs, 'PolynomA' takes priority over 'poly_a'.
        # - If it is present in kwargs, 'PolynomB' takes priority over 'poly_b' which
        #   in-turn takes priority over 'K' and 'H' if they are present in kwargs.
        # - Whenever this happens, we either raise an error or give a warning. If the
        #   duplicate definitions contain contradictory non-zero data then we raise an
        #   error, otherwise we give a warning.

        # Check kwargs and poly_a & poly_b for compatibility and convert to ParamArray
        prmpola = check_polynom("PolynomA", poly_a)
        prmpolb = check_polynom("PolynomB", poly_b)
        check_strength("K", 1)
        check_strength("H", 2)
        # Determine the length and order of PolynomA and PolynomB
        len_a, ord_a = getpol(prmpola)
        len_b, ord_b = getpol(prmpolb)
        deforder = max(getattr(self, "DefaultOrder", 0), ord_a, ord_b)
        # Remove MaxOrder
        maxorder = kwargs.pop("MaxOrder", deforder)
        kwargs.setdefault("PassMethod", "ThinMPolePass")
        super().__init__(family_name, **kwargs)
        # Set MaxOrder while PolynomA and PolynomB are not set yet
        super().__setattr__("MaxOrder", maxorder)
        # Adjust polynom lengths and set them
        len_ab = max(self.MaxOrder + 1, len_a, len_b)
        self.PolynomA = lengthen(prmpola, len_ab - len_a)
        self.PolynomB = lengthen(prmpolb, len_ab - len_b)

    def __setattr__(self, key, value):
        """Check the compatibility of MaxOrder, PolynomA and PolynomB"""
        polys = ("PolynomA", "PolynomB")
        if key in polys:
            lmin = self.MaxOrder
            if not len(value) > lmin:
                raise ValueError(f"Length of {key} must be larger than {lmin}")
        elif key == "MaxOrder":
            intval = int(value)
            lmax = min(len(getattr(self, k)) for k in polys)
            if not intval < lmax:
                raise ValueError(f"MaxOrder must be smaller than {lmax}")
        super().__setattr__(key, value)

    # noinspection PyPep8Naming
    @property
    def K(self) -> float:
        """Focusing strength [mˆ-2]"""
        arr = self.PolynomB
        return 0.0 if len(arr) < 2 else float(arr[1])

    # noinspection PyPep8Naming
    @K.setter
    def K(self, strength: float):
        self.PolynomB[1] = strength

    # noinspection PyPep8Naming
    @property
    def H(self) -> float:
        """Sextupolar strength [mˆ-3]"""
        arr = self.PolynomB
        return 0.0 if len(arr) < 3 else float(arr[2])

    # noinspection PyPep8Naming
    @H.setter
    def H(self, strength):
        self.PolynomB[2] = strength

    @property
    def IntegratedPolynomA(self) -> np.ndarray:
        """Integrated skew strengths."""
        return self.PolynomA

    @IntegratedPolynomA.setter
    def IntegratedPolynomA(self, value) -> None:
        self.PolynomA = value

    @property
    def IntegratedPolynomB(self) -> np.ndarray:
        """Integrated normal strengths."""
        return self.PolynomB

    @IntegratedPolynomB.setter
    def IntegratedPolynomB(self, value) -> None:
        self.PolynomB = value

    Kn0L = _k_property("IntegratedPolynomB", 0, "Integrated normal dipolar strength")
    Ks0L = _k_property("IntegratedPolynomA", 0, "Integrated skew dipolar strength")
    Kn1L = _k_property(
        "IntegratedPolynomB",
        1,
        "Integrated normal quadrupolar strength [m^-1]",
    )
    Ks1L = _k_property(
        "IntegratedPolynomA",
        1,
        "Integrated skew quadrupolar strength [m^-1]",
    )
    Kn2L = _k_property(
        "IntegratedPolynomB",
        2,
        "Integrated normal sextupolar strength in MAD convention [m^-2]",
    )
    Ks2L = _k_property(
        "IntegratedPolynomA",
        2,
        "Integrated skew sextupolar strength in MAD convention [m^-2]",
    )
    Kn3L = _k_property(
        "IntegratedPolynomB",
        3,
        "Integrated normal octupolar strength in MAD convention [m^-3]",
    )
    Ks3L = _k_property(
        "IntegratedPolynomA",
        3,
        "Integrated skew octupolar strength in MAD convention [m^-3]",
    )
    IntegratedStrength = _s_property(
        "IntegratedPolynomB",
        "Integrated strength of the main field component",
    )


class Multipole(_Radiative, LongElement, ThinMultipole):
    """Multipole element"""

    # Class attributes
    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["PolynomA", "PolynomB"]
    _stacklevel = 6  # Stacklevel for warnings

    # Instance attributes
    PolynomA: np.ndarray  #: Skew strengths
    PolynomB: np.ndarray  #: Straight strengths

    def __init__(self, family_name: str, length: float, poly_a, poly_b, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            poly_a:         Array of skew multipole components
            poly_b:         Array of normal multipole components

        Keyword arguments:
            MaxOrder:       Number of desired multipoles. Default: highest
              index of non-zero polynomial coefficients
            NumIntSteps:    Number of integration steps (default: 10)
            KickAngle:      Correction deviation angles (H, V)
            FieldScaling:   Scaling factor applied to the magnetic field
              (*PolynomA* and *PolynomB*)

        Default PassMethod: ``StrMPoleSymplectic4Pass``
        """
        kwargs.setdefault("PassMethod", "StrMPoleSymplectic4Pass")
        kwargs.setdefault("NumIntSteps", 10)
        super().__init__(family_name, length, poly_a, poly_b, **kwargs)

    def is_compatible(self, other) -> bool:
        if super().is_compatible(other) and self.MaxOrder == other.MaxOrder:
            for i in range(self.MaxOrder + 1):
                if self.PolynomB[i] != other.PolynomB[i]:
                    return False
                if self.PolynomA[i] != other.PolynomA[i]:
                    return False
            return True
        else:
            return False

    def _setter(self, attrname: str) -> Callable[[int, Any], None]:
        def setr(key, value):
            getattr(self, attrname)[key] = np.array(value) / self.Length

        return setr

    @property
    def IntegratedPolynomA(self) -> np.ndarray:
        """Integrated skew strengths (*Length* * *PolynomA*)."""
        pa = self.PolynomA
        ipa = _LinkedArray(pa.shape, dtype=pa.dtype, buffer=self.Length * pa)
        ipa.setitem = self._setter("PolynomA")
        return ipa

    @IntegratedPolynomA.setter
    def IntegratedPolynomA(self, value):
        self.PolynomA = np.array(value) / self.Length

    @property
    def IntegratedPolynomB(self) -> np.ndarray:
        """Integrated normal strengths (*Length* * *PolynomB*)."""
        pb = self.PolynomB
        ipb = _LinkedArray(pb.shape, dtype=pb.dtype, buffer=self.Length * pb)
        ipb.setitem = self._setter("PolynomB")
        return ipb

    @IntegratedPolynomB.setter
    def IntegratedPolynomB(self, value) -> None:
        self.PolynomB = np.array(value) / self.Length

    Strength = _s_property("PolynomB", "Strength of the main field component")
    Kn0 = _k_property("PolynomB", 0, "Normal dipolar strength [m^-1]")
    Ks0 = _k_property("PolynomA", 0, "Skew dipolar strength [m^-1]")
    Kn1 = _k_property("PolynomB", 1, "Normal quadrupolar strength [m^-2]")
    Ks1 = _k_property("PolynomA", 1, "Skew quadrupolar strength [m^-2]")
    Kn2 = _k_property(
        "PolynomB", 2, "Normal sextupolar strength in MAD convention [m^-3]"
    )
    Ks2 = _k_property(
        "PolynomA", 2, "Skew sextupolar strength in MAD convention [m^-3]"
    )
    Kn3 = _k_property(
        "PolynomB", 3, "Normal octupolar strength in MAD convention [m^-4]"
    )
    Ks3 = _k_property("PolynomA", 3, "Skew octupolar strength in MAD convention [m^-4]")


class Dipole(Radiative, Multipole):
    """Dipole element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["BendingAngle", "K"]
    _conversions = dict(
        Multipole._conversions,
        EntranceAngle=float,
        ExitAngle=float,
        FringeInt1=float,
        FringeInt2=float,
        FringeQuadEntrance=int,
        FringeQuadExit=int,
        FringeBendEntrance=int,
        FringeBendExit=int,
    )
    _stacklevel = 7  # Stacklevel for warnings

    _entrance_fields = Multipole._entrance_fields + [
        "EntranceAngle",
        "FringeInt1",
        "FringeBendEntrance",
        "FringeQuadEntrance",
    ]
    _exit_fields = Multipole._exit_fields + [
        "ExitAngle",
        "FringeInt2",
        "FringeBendExit",
        "FringeQuadExit",
    ]

    DefaultOrder = 0

    def __init__(
        self,
        family_name: str,
        length: float,
        bending_angle: float = 0.0,
        k: float = 0.0,
        **kwargs,
    ):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            bending_angle:  Bending angle [rd]
            k:              Focusing strength [m^-2]

        Keyword arguments:
            EntranceAngle=0.0:  entrance angle
            ExitAngle=0.0:      exit angle
            PolynomB:           straight multipoles
            PolynomA:           skew multipoles
            MaxOrder=0:         Number of desired multipoles
            NumIntSt=10:        Number of integration steps
            FullGap:            Magnet full gap
            FringeInt1:         Extension of the entrance fringe field
            FringeInt2:         Extension of the exit fringe field
            FringeBendEntrance: 1: legacy version Brown First Order (default)

              2: SOLEIL close to second order of Brown

              3: THOMX
            FringeBendExit:     See *FringeBendEntrance*
            FringeQuadEntrance: 0: no fringe field effect (default)

              1: Lee-Whiting's thin lens limit formula

              2: elegant-like
            FringeQuadExit:     See *FringeQuadEntrance*
            fringeIntM0:        Integrals for FringeQuad method 2
            fringeIntP0:
            KickAngle:          Correction deviation angles (H, V)
            FieldScaling:       Scaling factor applied to the magnetic field

        Available PassMethods: :ref:`BndMPoleSymplectic4Pass`,
        :ref:`BendLinearPass`, :ref:`ExactSectorBendPass`,
        :ref:`ExactRectangularBendPass`, :ref:`ExactRectBendPass`,
        BndStrMPoleSymplectic4Pass

        Default PassMethod: :ref:`BndMPoleSymplectic4Pass`
        """
        kwargs.setdefault("BendingAngle", bending_angle)
        kwargs.setdefault("EntranceAngle", 0.0)
        kwargs.setdefault("ExitAngle", 0.0)
        kwargs.setdefault("PassMethod", "BndMPoleSymplectic4Pass")
        super().__init__(family_name, length, [], [0.0, k], **kwargs)

    def items(self) -> Generator[tuple[str, Any], None, None]:
        yield from super().items()
        yield "K", self.K

    def _part(self, fr, sumfr):
        pp = super()._part(fr, sumfr)
        pp.BendingAngle = fr / sumfr * self.BendingAngle
        pp.EntranceAngle = 0.0
        pp.ExitAngle = 0.0
        return pp

    def is_compatible(self, other) -> bool:
        def invrho(dip: Dipole):
            return dip.BendingAngle / dip.Length

        return (
            super().is_compatible(other)
            and self.ExitAngle == -other.EntranceAngle
            and abs(invrho(self) - invrho(other)) <= 1.0e-6
        )

    def merge(self, other) -> None:
        super().merge(other)
        # noinspection PyAttributeOutsideInit
        self.ExitAngle = other.ExitAngle
        self.BendingAngle += other.BendingAngle


# Bend is a synonym of Dipole.
Bend = Dipole


class Quadrupole(Radiative, Multipole):
    """Quadrupole element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["K"]
    _conversions = dict(
        Multipole._conversions, FringeQuadEntrance=int, FringeQuadExit=int
    )
    _stacklevel = 7  # Stacklevel for warnings

    _entrance_fields = Multipole._entrance_fields + ["FringeQuadEntrance"]
    _exit_fields = Multipole._exit_fields + ["FringeQuadExit"]

    DefaultOrder = 1

    def __init__(self, family_name: str, length: float, k: float = 0.0, **kwargs):
        """Quadrupole(FamName, Length, Strength=0, **keywords)

        Args:
            family_name:    Name of the element
            length:         Element length [m]
            k:              Focusing strength [mˆ-2]

        Keyword Arguments:
            PolynomB:           straight multipoles
            PolynomA:           skew multipoles
            MaxOrder=1:         Number of desired multipoles
            NumIntSteps=10:     Number of integration steps
            FringeQuadEntrance: 0: no fringe field effect (default)

              1: Lee-Whiting's thin lens limit formula

              2: elegant-like
            FringeQuadExit:     See ``FringeQuadEntrance``
            fringeIntM0:        Integrals for FringeQuad method 2
            fringeIntP0:
            KickAngle:          Correction deviation angles (H, V)
            FieldScaling:       Scaling factor applied to the magnetic field
              (*PolynomA* and *PolynomB*)

        Default PassMethod: ``StrMPoleSymplectic4Pass``
        """
        kwargs.setdefault("PassMethod", "StrMPoleSymplectic4Pass")
        super().__init__(family_name, length, [], [0.0, k], **kwargs)

    def items(self) -> Generator[tuple[str, Any], None, None]:
        yield from super().items()
        yield "K", self.K


class Sextupole(Multipole):
    """Sextupole element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["H"]
    _stacklevel = 7  # Stacklevel for warnings

    DefaultOrder = 2

    def __init__(self, family_name: str, length: float, h: float = 0.0, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            h:              Sextupolar strength [mˆ-3]

        Keyword Arguments:
            PolynomB:           straight multipoles
            PolynomA:           skew multipoles
            MaxOrder:           Number of desired multipoles
            NumIntSteps=10:     Number of integration steps
            KickAngle:          Correction deviation angles (H, V)
            FieldScaling:       Scaling factor applied to the magnetic field
              (*PolynomA* and *PolynomB*)

        Default PassMethod: ``StrMPoleSymplectic4Pass``
        """
        kwargs.setdefault("PassMethod", "StrMPoleSymplectic4Pass")
        super().__init__(family_name, length, [], [0.0, 0.0, h], **kwargs)


class Octupole(Multipole):
    """Octupole element, with no changes from multipole at present"""

    _BUILD_ATTRIBUTES = Multipole._BUILD_ATTRIBUTES

    DefaultOrder = 3


class Corrector(LongElement):
    """Corrector element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["KickAngle"]

    def __init__(self, family_name: str, length: float, kick_angle, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            KickAngle:      Correction deviation angles (H, V)

        Keyword Args:
            FieldScaling:   Scaling factor applied to the magnetic field
              (*KickAngle*)

        Default PassMethod: ``CorrectorPass``
        """
        kwargs.setdefault("PassMethod", "CorrectorPass")
        super().__init__(family_name, length, KickAngle=kick_angle, **kwargs)


class Wiggler(Radiative, LongElement):
    """Wiggler element

    See atwiggler.m
    """

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["Lw", "Bmax"]
    _conversions = dict(
        Element._conversions,
        Lw=_float,
        Bmax=_float,
        Energy=_float,
        Bx=lambda v: _array(v, (6, -1)),
        By=lambda v: _array(v, (6, -1)),
        Nstep=int,
        Nmeth=int,
        NHharm=int,
        NVharm=int,
    )

    # noinspection PyPep8Naming
    def __init__(
        self,
        family_name: str,
        length: float,
        wiggle_period: float,
        b_max: float,
        energy: float = 0.0,
        *,
        Nstep: int | None = 5,
        Nmeth: int | None = 4,
        By=(1, 1, 0, 1, 1, 0),
        Bx=(),
        **kwargs,
    ):
        """
        Args:
            length:         total length of the wiggler
            wiggle_period:  length must be a multiple of this
            b_max:          peak wiggler field [Tesla]
            energy:         kept for backwards compatibility but always ignored
            Nstep:          number of integration steps.
            Nmeth:          symplectic integration order: 2 or 4
            Bx:             harmonics for horizontal wiggler: (6, nHharm)
                              array-like object
            By:             harmonics for vertical wiggler (6, nHharm)
                              array-like object

        Default PassMethod: ``GWigSymplecticPass``
        """
        kwargs.setdefault("PassMethod", "GWigSymplecticPass")
        n_wiggles = length / wiggle_period
        if abs(round(n_wiggles) - n_wiggles) > 1e-6:
            raise ValueError(
                "Wiggler: length / wiggle_period is not an "
                f"integer. ({length}/{wiggle_period}={n_wiggles})"
            )
        if energy:
            warnings.warn(
                AtWarning(
                    f"Wiggler: positional argument energy={energy} has been ignored as "
                    "lattice.Energy is always used instead."
                )
            )
        super().__init__(
            family_name,
            length,
            Lw=wiggle_period,
            Bmax=b_max,
            Nstep=Nstep,
            Nmeth=Nmeth,
            By=By,
            Bx=Bx,
            **kwargs,
        )

        for i, b in enumerate(self.By.T):
            dk = abs(b[3] ** 2 - b[4] ** 2 - b[2] ** 2) / abs(b[4])
            if dk > 1e-6:
                raise ValueError(f"Wiggler(H): kx^2 + kz^2 -ky^2 !=0, i = {i}")

        for i, b in enumerate(self.Bx.T):
            dk = abs(b[2] ** 2 - b[4] ** 2 - b[3] ** 2) / abs(b[4])
            if dk > 1e-6:
                raise ValueError(f"Wiggler(V): ky^2 + kz^2 -kx^2 !=0, i = {i}")

        self.NHharm = self.By.shape[1]
        self.NVharm = self.Bx.shape[1]

    def divide(self, frac) -> list[Element]:
        # A wiggler is indivisible
        return [self]
