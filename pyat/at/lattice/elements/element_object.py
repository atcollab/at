"""Base :py:class:`.Element` object"""

from __future__ import annotations

__all__ = ["Element"]

import re
from collections.abc import Generator
from copy import copy, deepcopy
from typing import Any

import numpy as np

from .conversions import _array, _array66, _int, _float
from ..parambase import ParamDef, _nop
from ..parameters import _ACCEPTED, Param, ParamArray


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
            else:
                value = conversion(value)
        except Exception as exc:
            # Conversion failed
            exc.args = (f"{self._ident(attrname)}: {exc}",)
            raise
        else:
            # Set the attribute with the converted value
            object.__setattr__(self, attrname, value)

    def __getattribute__(self, attrname: str) -> Any:
        """Override __getattribute__ to handle parameter values.

        This method returns the value of parameters instead of the parameter objects
        themselves when accessing attributes.
        """
        try:
            attr = object.__getattribute__(self, attrname)
        except AttributeError as exc:
            cl = self.__class__.__name__
            el = object.__getattribute__(self, "FamName")
            exc.args = (f"{cl}({el!r}) has no attribute {attrname!r}",)
            raise

        # If it's a parameter or parameter array, return its value
        if isinstance(attr, (ParamDef, ParamArray)):
            return attr.value

        # Otherwise return the attribute itself
        return attr

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
        return set(vars(self).keys())

    def to_dict(self, freeze: bool = True):
        """Return a copy of the element parameters"""
        if freeze:
            return {k: getattr(self, k) for k in self.keys()}
        else:
            return vars(self).copy()

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
        """split the element in len(frac) pieces whose length is frac[i]*self.Length

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

    def set_parameter(
        self, attrname: str, value: Any, index: int | None = None
    ) -> None:
        """Set an element's parameter.

        This allows setting a parameter into an attribute or an item of an
        array attribute.

        Args:
            attrname:   Attribute name
            value:      Parameter or value to be set
            index:      Index into an array attribute. If *value* is a
              parameter, the array attribute is converted to a :py:class:`.ParamArray`.

        Raises:
            IndexError: If the provided index is out of bounds for the array attribute
            AttributeError: If the attribute doesn't exist
        """

        def set_array_item(arr: np.ndarray, idx: int, val: Any) -> None:
            """Helper function to set an item in an array."""
            try:
                arr[idx] = val
            except IndexError as exc:
                exc.args = (f"{self._ident(attrname)}: {exc}",)
                raise

        # Set the entire attribute
        if index is None:
            setattr(self, attrname, value)
        # Set a specific index in an array attribute
        else:
            array = self.get_parameter(attrname)
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
        # Check AT or MADX parameters
        # Check if any attribute is parameterised
        if attrname is None:
            return any(self.is_parameterised(attribute) for attribute in self.__dict__)

        # Get the attribute or specific index
        attribute = self.get_parameter(attrname, index=index)

        # Check if the attribute itself is a parameter
        if isinstance(attribute, ParamDef):
            return True

        # Check if any element in the array is a parameter
        if isinstance(attribute, np.ndarray):
            return any(isinstance(item, ParamDef) for item in attribute.flat)

        # Not parameterised
        return False

    def parameterise(
        self, attrname: str, index: int | None = None, name: str = ""
    ) -> _ACCEPTED:
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
        if isinstance(current_value, _ACCEPTED):
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

    def unparameterise(
        self, attrname: str | None = None, index: int | None = None
    ) -> None:
        """Replace parameters with their current values.

        This function replaces parameters with their current values, effectively
        "freezing" them. This is useful when you want to convert a parameterised
        element back to a regular element with fixed values.

        Args:
            attrname:   Attribute name. If :py:obj:`None`, freezes all attributes
            index:      Index in an array. If :py:obj:`None`, freezes the whole
              attribute

        Attributes which are not parameters are silently ignored.
        """

        def _freeze_attribute(attrname: str, attr: Any) -> None:
            """Helper function to freeze a parameterised attribute."""
            # Accepts AT or MADX parameters
            if isinstance(attr, ParamDef):
                setattr(self, attrname, attr.value)
            elif isinstance(attr, np.ndarray):
                for i, item in enumerate(attr.flat):
                    if isinstance(item, ParamDef):
                        ij = np.unravel_index(i, attr.shape)
                        attr[ij] = item.value

        if attrname is None:
            for name, attr in self.__dict__.items():
                _freeze_attribute(name, attr)
        else:
            attr = self.get_parameter(attrname)
            if index is None:
                _freeze_attribute(attrname, attr)
            else:
                item = attr[index]
                if isinstance(item, ParamDef):
                    attr[index] = item.value
