from __future__ import annotations

__all__ = [
    "set_parameter",
    "is_parameterised",
    "parameterise",
    "unparameterise",
]

from typing import Any

import numpy as np

from .elements import Element
from .parser import ParamDef, _nop
from .variables import ParamBase
from .parameters import Param, ParamArray, _ACCEPTED


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
            # Convert the array to a ParamArray if it's not already one'
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


def _setattr(self, attrname: str, value: Any) -> None:
    """Override __setattr__ to handle parameter conversions.

    This method applies the appropriate conversion function to the value
    before setting it as an attribute.

    Args:
        attrname: Name of the attribute to set
        value: Value to set for the attribute

    Raises:
        Exception: If the conversion fails
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


def _getattr(self, attrname) -> Any:
    """Override __getattr__ to handle parameter values.

    This method returns the value of parameters instead of the parameter objects
    themselves when accessing attributes.

    Args:
        attrname: Name of the attribute to get

    Returns:
        The attribute value, or the parameter value if the attribute is a parameter
    """
    try:
        return self._parameters[attrname].value
    except KeyError as exc:
        cl = self.__class__.__name__
        el = object.__getattribute__(self, "FamName")
        raise AttributeError(f"{cl}({el!r}) has no attribute {attrname!r}") from exc


def _delattr(self, attrname) -> None:
    try:
        del self._parameters[attrname]
    except KeyError:
        object.__delattr__(self, attrname)


Element.__setattr__ = _setattr
Element.__getattr__ = _getattr
Element.__delattr__ = _delattr
Element.set_parameter = set_parameter
Element.is_parameterised = is_parameterised
Element.parameterise = parameterise
Element.unparameterise = unparameterise
