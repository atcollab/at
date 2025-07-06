"""
Module to define common elements used in AT.

Each element has a default PassMethod attribute for which it should have the
appropriate attributes.  If a different PassMethod is set, it is the caller's
responsibility to ensure that the appropriate attributes are present.
"""

from __future__ import annotations

__all__ = [
    "LongtMotion",
    "_DictLongtMotion",
    "_Radiative",
    "Radiative",
    "Collective",
    "Element",
    "_array",
    "_array66",
    "_float",
    "_int",
]

import abc
import re
import warnings
from abc import ABC
from collections.abc import Generator
from copy import copy, deepcopy
from typing import Any

import numpy as np

from .exceptions import AtWarning
from .parser import _nop, ParamDef
from .variables import ParamBase
from .parameters import _ACCEPTED, Param, ParamArray, AttributeArray as _array

# AtWarning from this module should always be issued (not only on the first occurrence)
warnings.filterwarnings("always", category=AtWarning, module=__name__)


def _array66(value):
    return _array(value, shape=(6, 6))


def _float(value) -> float:
    return float(value)


def _int(value, vmin: int | None = None, vmax: int | None = None) -> int:
    intv = int(value)
    if vmin is not None and intv < vmin:
        raise ValueError(f"Value must be greater of equal to {vmin}")
    if vmax is not None and intv > vmax:
        raise ValueError(f"Value must be smaller of equal to {vmax}")
    return intv


class LongtMotion(ABC):
    """Abstract Base class for all Element classes whose instances may modify
    the particle momentum

    Allows identifying elements potentially inducing longitudinal motion.

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
            enable:     :py:obj:`True`: for enabling, :py:obj:`False` for
              disabling
            new_pass:   New PassMethod:

              * :py:obj:`None`: makes no change,
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
        return None


# noinspection PyUnresolvedReferences
class _DictLongtMotion(LongtMotion):
    # noinspection PyShadowingNames
    """Mixin class for elements implementing a 'default_pass' class attribute

    :py:class:`DictLongtMotion` provides:

    * a :py:meth:`set_longt_motion` method setting the PassMethod according
      to the ``default_pass`` dictionary.
    * a :py:obj:`.longt_motion` property set to :py:obj:`True` when the
      PassMethod is ``default_pass[True]``

    The class must have a ``default_pass`` class attribute, a dictionary
    such that:

    * ``default_pass[False]`` is the PassMethod when radiation is turned
      OFF,
    * ``default_pass[True]`` is the default PassMethod when radiation is
      turned ON.

    The :py:class:`DictLongtMotion` class must be set as the first base class.

    Example:

        >>> class QuantumDiffusion(_DictLongtMotion, Element):
        ...     default_pass = {False: "IdentityPass", True: "QuantDiffPass"}

        Defines a class such that :py:meth:`set_longt_motion` will select
        ``'IdentityPass'`` or ``'IdentityPass'``.
    """

    def _get_longt_motion(self):
        return self.PassMethod != self.default_pass[False]

    # noinspection PyShadowingNames
    def set_longt_motion(self, enable, new_pass=None, **kwargs):
        if new_pass == "auto":
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
        return self.PassMethod.endswith(("RadPass", "QuantPass"))

    def _autopass(self, enable):
        if enable:
            root = self.PassMethod.replace("QuantPass", "Pass").replace(
                "RadPass", "Pass"
            )
            return "".join((root[:-4], "RadPass"))
        elif self.longt_motion:
            root = self.PassMethod.replace("QuantPass", "Pass").replace(
                "RadPass", "Pass"
            )
            return root
        else:
            return None

    # noinspection PyTypeChecker,PyShadowingNames
    def set_longt_motion(self, enable, new_pass=None, copy=False, **kwargs):
        if new_pass == "auto":
            new_pass = self._autopass(enable)
        if new_pass is None or new_pass == self.PassMethod:
            return self if copy else None
        if enable:

            def setpass(el):
                el.PassMethod = new_pass
                el.Energy = kwargs["energy"]

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
        return None


class Radiative(_Radiative):
    # noinspection PyUnresolvedReferences
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

    Derived classes will automatically set the
    :py:attr:`~Element.is_collective` property when the element is active.

    The class must have a ``default_pass`` class attribute, a dictionary such
    that:

    * ``default_pass[False]`` is the PassMethod when collective effects
      are turned OFF,
    * ``default_pass[True]`` is the default PassMethod when collective effects
      are turned ON.

    The :py:class:`Collective` class must be set as the first base class.

    Example:
        >>> class WakeElement(Collective, Element):
        ...     default_pass = {False: "IdentityPass", True: "WakeFieldPass"}

        Defines a class where the :py:attr:`~Element.is_collective` property is
        handled
    """

    def _get_collective(self):
        # noinspection PyUnresolvedReferences
        return self.PassMethod != self.default_pass[False]

    @abc.abstractmethod
    def clear_history(self):
        pass


class Element:
    """Base class for AT elements"""

    _BUILD_ATTRIBUTES = ["FamName"]
    _conversions = {
        "FamName": str,
        "PassMethod": str,
        "Length": _float,
        "R1": _array66,
        "R2": _array66,
        "T1": lambda v: _array(v, (6,)),
        "T2": lambda v: _array(v, (6,)),
        "RApertures": lambda v: _array(v, (4,)),
        "EApertures": lambda v: _array(v, (2,)),
        "KickAngle": lambda v: _array(v, (2,)),
        "PolynomB": _array,
        "PolynomA": _array,
        "BendingAngle": _float,
        "MaxOrder": _int,
        "NumIntSteps": lambda v: _int(v, vmin=0),
        "Energy": _float,
    }

    _entrance_fields = ["T1", "R1"]
    _exit_fields = ["T2", "R2"]
    _no_swap = _entrance_fields + _exit_fields
    __slots__ = ["_parameters", "__dict__"]

    def __new__(cls, *args, **kwargs):
        obj = super().__new__(cls)
        # _parameters must be created before any other attribute is set
        obj._parameters = {}
        return obj

    def __init__(self, family_name: str, **kwargs):
        """
        Parameters:
            family_name:    Name of the element

        All keywords will be set as attributes of the element
        """

        self.FamName = family_name
        self.Length = kwargs.pop("Length", 0.0)
        self.PassMethod = kwargs.pop("PassMethod", "IdentityPass")
        self.update(kwargs)

    def __setattr__(self, attrname: str, value: Any) -> None:
        """Override __setattr__ to handle parameter conversions.

        This method applies the appropriate conversion function to the value
        before setting it as an attribute.
        """
        # Get the conversion function for this attribute or use _nop (no operation)
        conversion = self._conversions.get(attrname, _nop)

        try:
            # If the value is a parameter, set its conversion function
            if isinstance(value, _ACCEPTED):
                value.set_conversion(conversion)
            # Otherwise, apply the conversion to the value
            elif not isinstance(value, ParamArray):
                value = conversion(value)
        except Exception as exc:
            # Conversion failed
            exc.args = (f"{self._ident(attrname)}: {exc}",)
            raise
        else:
            # Conversion succeeded
            if isinstance(value, (ParamDef, ParamArray)):
                # Store the parameter and remove the attribute
                self._parameters[attrname] = value
                try:
                    object.__delattr__(self, attrname)
                except AttributeError:
                    # Attribute doesn't exist, which is fine
                    pass
            else:
                # Store the attribute and remove the parameter
                object.__setattr__(self, attrname, value)
                try:
                    del self._parameters[attrname]
                except KeyError:
                    # Parameter doesn't exist, which is fine
                    pass

    def __getattr__(self, attrname: str) -> Any:
        """Override __getattr__ to handle parameter values.

        This method returns the value of parameters instead of the parameter objects
        themselves when accessing attributes.
        """
        try:
            return self._parameters[attrname].value
        except KeyError as exc:
            cl = self.__class__.__name__
            el = object.__getattribute__(self, "FamName")
            raise AttributeError(f"{cl}({el!r}) has no attribute {attrname!r}") from exc

    def __delattr__(self, attrname: str) -> None:
        """Override __delattr__ to handle parameter deletions."""
        try:
            object.__delattr__(self, attrname)
        except AttributeError:
            try:
                del self._parameters[attrname]
            except KeyError as exc:
                cl = self.__class__.__name__
                el = object.__getattribute__(self, "FamName")
                raise AttributeError(
                    f"{cl}({el!r}) has no attribute {attrname!r}"
                ) from exc

    def __str__(self):
        return "\n".join(
            [self.__class__.__name__ + ":"]
            + [f"{k:>14}: {v!s}" for k, v in self.items(freeze=False)]
        )

    def __repr__(self):
        clsname, args, kwargs = self.definition
        keywords = [f"{arg!r}" for arg in args]
        keywords += [f"{k}={v!r}" for k, v in kwargs.items()]
        args = re.sub(r"\n\s*", " ", ", ".join(keywords))
        return f"{clsname}({args})"

    def __getstate__(self):
        # Make a copy of parameters
        return self.__dict__, {"_parameters": self._parameters.copy()}

    def _ident(self, attrname: str | None = None, index: bool = None):
        """Return an element's identifier for error messages"""
        if attrname is None:
            return f"{self.__class__.__name__}({self.FamName!r})"
        elif index is None:
            return f"{self.__class__.__name__}({self.FamName!r}).{attrname}"
        else:
            return f"{self.__class__.__name__}({self.FamName!r}).{attrname}[{index}]"

    @classmethod
    def subclasses(cls) -> Generator[type[Element], None, None]:
        """Yields all the class subclasses.

        Some classes may appear several times because of diamond-shape inheritance
        """
        for subclass in cls.__subclasses__():
            yield from subclass.subclasses()
        yield cls

    def keys(self):
        """Return a set of all attribute names"""
        v = set(vars(self).keys())
        v.update(self._parameters.keys())
        return v

    def to_dict(self, freeze: bool = True):
        """Return a copy of the element parameters"""
        if freeze:
            return {k: getattr(self, k) for k in self.keys()}
        else:
            v = vars(self).copy()
            v.update(self._parameters)
            return v

    def get_parameter(self, attrname: str, index: int | None = None) -> Any:
        """Extract a parameter of an element

        Unlike :py:func:`getattr`, :py:func:`get_parameter` returns the
        parameter itself instead of its value. If the item is not a parameter,
        both functions are equivalent, the value is returned.

        Args:
            attrname:   Attribute name
            index:      Index in an array attribute. If :py:obj:`None`, the
              whole attribute is returned

        Returns:
            The parameter object or attribute value.
        """
        try:
            attr = self.__dict__[attrname]
        except KeyError:
            try:
                attr = self._parameters[attrname]
            except KeyError:
                raise AttributeError(
                    f"{self._ident()} has no attribute {attrname!r}"
                ) from None
        if index is not None:
            try:
                attr = attr[index]
            except IndexError as exc:
                raise IndexError(f"{self._ident(attrname)}: {exc}") from None
        return attr

    def equals(self, other) -> bool:
        """Whether an element is equivalent to another.

        This implementation was found to be too slow for the generic
        __eq__ method when comparing lattices.
        """
        return repr(self) == repr(other)

    def divide(self, frac) -> list[Element]:
        # noinspection PyUnresolvedReferences
        """split the element in len(frac) pieces whose length
        is frac[i]*self.Length

        Parameters:
            frac:           length of each slice expressed as a fraction of the
              initial length. ``sum(frac)`` may differ from 1.

        Returns:
            elem_list:  a list of elements equivalent to the original.

        Example:

            >>> Drift("dr", 0.5).divide([0.2, 0.6, 0.2])
            [Drift('dr', 0.1), Drift('dr', 0.3), Drift('dr', 0.1)]
        """
        # Bx default, the element is indivisible
        return [self]

    def swap_faces(self, copy=False):
        """Swap the faces of an element, alignment errors are ignored"""

        def swapattr(element, attro, attri):
            val = element.get_parameter(attri)  # get the parameter itself
            delattr(element, attri)
            return attro, val

        if copy:
            el = self.copy()
        else:
            el = self
        # Remove and swap entrance and exit attributes
        attrs = el.keys()
        fin = dict(
            swapattr(el, kout, kin)
            for kin, kout in zip(el._entrance_fields, el._exit_fields)
            if kin in attrs and kin not in el._no_swap
        )
        fout = dict(
            swapattr(el, kin, kout)
            for kin, kout in zip(el._entrance_fields, el._exit_fields)
            if kout in attrs and kout not in el._no_swap
        )
        # Apply swapped entrance and exit attributes
        for key, value in fin.items():
            setattr(el, key, value)
        for key, value in fout.items():
            setattr(el, key, value)
        return el if copy else None

    def update(self, *args, **kwargs):
        """
        update(**kwargs)
        update(mapping, **kwargs)
        update(iterable, **kwargs)
        Update the element attributes with the given arguments
        """
        attrs = dict(*args, **kwargs)
        for key, value in attrs.items():
            setattr(self, key, value)

    def copy(self) -> Element:
        """Return a shallow copy of the element"""
        return copy(self)

    def deepcopy(self) -> Element:
        """Return a deep copy of the element"""
        return deepcopy(self)

    @property
    def definition(self) -> tuple[str, tuple, dict]:
        """tuple (class_name, args, kwargs) defining the element"""
        attrs = {k: getattr(self, k) for k, v in self.items()}
        arguments = tuple(
            attrs.pop(k, getattr(self, k)) for k in self._BUILD_ATTRIBUTES
        )
        defelem = self.__class__(*arguments)
        keywords = {
            k: v
            for k, v in attrs.items()
            if not np.array_equal(v, getattr(defelem, k, None))
        }
        return self.__class__.__name__, arguments, keywords

    def items(self, freeze: bool = True) -> Generator[tuple[str, Any], None, None]:
        """Iterates through the data members"""
        v = self.to_dict(freeze=freeze)
        for k in ["FamName", "Length", "PassMethod"]:
            yield k, v.pop(k)
        for k, val in sorted(v.items()):
            yield k, val

    def is_compatible(self, other: Element) -> bool:
        """Checks if another :py:class:`Element` can be merged"""
        return False

    def merge(self, other) -> None:
        """Merge another element"""
        if not self.is_compatible(other):
            badname = getattr(other, "FamName", type(other))
            raise TypeError(f"Cannot merge {self.FamName} and {badname}")

    # noinspection PyMethodMayBeStatic
    def _get_longt_motion(self):
        return False

    # noinspection PyMethodMayBeStatic
    def _get_collective(self):
        return False

    @property
    def longt_motion(self) -> bool:
        """:py:obj:`True` if the element affects the longitudinal motion"""
        return self._get_longt_motion()

    @property
    def is_collective(self) -> bool:
        """:py:obj:`True` if the element involves collective effects"""
        return self._get_collective()

    def set_parameter(self, attrname: str, value: Any, index: int | None = None) -> None:
        """Set an element's parameter.

        This allows setting a parameter into an attribute or an item of an array attribute.

        Args:
            attrname:   Attribute name
            value:      Parameter or value to be set
            index:      Index into an array attribute. If *value* is a
              parameter, the array attribute is converted to a :py:class:`.ParamArray`.

        Raises:
            IndexError: If the provided index is out of bounds for the array attribute
            AttributeError: If the attribute doesn't exist
        """

        def set_array_item(array: np.ndarray, idx: int, val: Any) -> None:
            """Helper function to set an item in an array with improved error handling."""
            try:
                array[idx] = val
            except IndexError as exc:
                exc.args = (f"{self._ident(attrname)}: {exc}",)
                raise

        # Set the entire attribute
        if index is None:
            setattr(self, attrname, value)
        # Set a specific index in an array attribute
        else:
            array = self.get_parameter(attrname)
            if not isinstance(array, ParamArray) and isinstance(value, ParamDef):
                # Convert the array to a ParamArray if it's not already one
                array = ParamArray(array, shape=array.shape, dtype=array.dtype)
                set_array_item(array, index, value)
                setattr(self, attrname, array)
            else:
                set_array_item(array, index, value)


    def is_parameterised(
        self, attrname: str | None = None, index: int | None = None
    ) -> bool:
        """Check for the parameterisation of an element

        Args:
            attrname:   Attribute name. If :py:obj:`None`, checks if any attribute is
              parameterised
            index:      Index in an array attribute. If :py:obj:`None`, tests the whole
              attribute

        Returns:
            True if the attribute, or array item is parameterised, False otherwise
        """
        # Check if any attribute is parameterised
        if attrname is None:
            return len(self._parameters) > 0

        # Get the attribute or specific index
        attribute = self.get_parameter(attrname, index=index)
        return isinstance(attribute, (ParamDef, ParamArray))


    def parameterise(
        self, attrname: str, index: int | None = None, name: str = ""
    ) -> ParamBase:
        """Convert attribute to parameter preserving value.

        The value of the attribute is kept unchanged. If the attribute is
        already parameterised, the existing parameter is returned.

        Args:
            attrname:   Attribute name
            index:      Index in an array. If :py:obj:`None`, the
              whole attribute is parameterised
            name:       Name of the created parameter

        Returns:
            A :py:class:`.ParamArray` for an array attribute,
            a :py:class:`.Param` for a scalar attribute or an item in an
            array attribute

        Raises:
            TypeError: If the attribute value is not a valid parameter type (Number)
            IndexError: If the index is out of bounds for an array attribute
            AttributeError: If the attribute does not exist
        """
        # Get the current value of the attribute or array element
        current_value = self.get_parameter(attrname, index=index)

        # If it's already a parameter, return it
        if isinstance(current_value, ParamBase):
            return current_value

        # Create a new parameter with the current value
        try:
            param = Param(current_value, name=name)
        except TypeError as exc:
            exc.args = (f"Cannot parameterise {self._ident(attrname)}: {exc}",)
            raise

        # Set the parameter in the element
        self.set_parameter(attrname, param, index=index)
        return param


    def unparameterise(self, attrname: str | None = None, index: int | None = None) -> None:
        """Freeze the parameter values by replacing parameters with their current values.

        This function replaces parameters with their current values, effectively
        "freezing" them. This is useful when you want to convert a parameterised
        element back to a regular element with fixed values.

        Args:
            attrname:   Attribute name. If :py:obj:`None`, freezes all attributes
            index:      Index in an array. If :py:obj:`None`, freezes the whole attribute

        Attributes which are not parameters are silently ignored.
        """
        if attrname is None:
            # freeze all the attributes
            # make a copy of the parameters dict to avoid modifications during iteration
            for name, attr in self._parameters.copy().items():
                setattr(self, name, attr.value)
        else:
            attr = self.get_parameter(attrname)
            if not isinstance(attr, (ParamDef, ParamArray)):
                # silently ignore non-parameter attributes
                return
            if index is None:
                # freeze a scalar attribute
                setattr(self, attrname, attr.value)
            else:
                # freeze an item in an array attribute
                item = attr[index]
                if isinstance(item, ParamDef):
                    attr[index] = item.value
                if not any(isinstance(item, ParamDef) for item in attr.flat):
                    # freeze the whole array attribute if none of its items is a parameter
                    setattr(self, attrname, attr.value)
