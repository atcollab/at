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


def set_parameter(self, attrname: str, value, index: int | None = None) -> None:
    """Set an element's parameter

    This allows setting a parameter into an item of an array attribute.

    Args:
        attrname:   Attribute name
        value:      Parameter or value to be set
        index:      Index into an array attribute. If *value* is a
          parameter, the array attribute is converted to a :py:class:`.ParamArray`.
    """

    def setitem(array, index, value):
        try:
            array[index] = value
        except IndexError as exc:
            exc.args = (f"{self._ident(attrname)}: {exc}",)
            raise

    if index is None:
        setattr(self, attrname, value)
    else:
        array = self.get_parameter(attrname)
        setitem(array, index, value)


def is_parameterised(
    self, attrname: str | None = None, index: int | None = None
) -> bool:
    """Check for the parametrisation of an element

    Args:
        attrname:   Attribute name. If :py:obj:`None`, return :py:obj:`True`
          if any attribute is parameterised
        index:      Index in an array attribute. If :py:obj:`None`, the
          whole attribute is tested for parameterisation

    Returns:
        True if the attribute, or array item is parameterised, False otherwise
    """
    if attrname is None:
        for attr in self.__dict__:
            if self.is_parameterised(attr):
                return True
        return False
    else:
        attr = self.get_parameter(attrname, index=index)
        if isinstance(attr, ParamDef):
            return True
        elif isinstance(attr, np.ndarray):
            return any(isinstance(item, ParamDef) for item in attr.flat)
        else:
            return False


def parameterise(
    self, attrname: str, index: int | None = None, name: str = ""
) -> ParamBase:
    """Convert an attribute into a parameter

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
    vini = self.get_parameter(attrname, index=index)

    if isinstance(vini, ParamBase):
        return vini

    try:
        param = Param(vini, name=name)
    except TypeError as exc:
        exc.args = (f"Cannot parameterise {self._ident(attrname)}: {exc}",)
        raise

    self.set_parameter(attrname, param, index=index)
    return param


def unparameterise(self, attrname: str | None = None, index: int | None = None) -> None:
    """Freeze the parameter values by replacing parameters with their current values.

    This function replaces parameters with their current values, effectively
    "freezing" them. This is useful when you want to convert a parameterised
    element back to a regular element with fixed values.

    Args:
        attrname:   Attribute name. If :py:obj:`None`, all the attributes
          are frozen
        index:      Index in an array. If :py:obj:`None`, the whole
          attribute is frozen

    Attributes which are not parameters are silently ignored.
    """

    def unparam_attr(attrname: str, attr: Any) -> None:
        """Helper function to unparameterise a single attribute."""
        if isinstance(attr, ParamDef):
            setattr(self, attrname, attr.value)
        elif isinstance(attr, np.ndarray):
            for i, item in enumerate(attr.flat):
                if isinstance(item, ParamDef):
                    ij = np.unravel_index(i, attr.shape)
                    attr[ij] = item.value

    if attrname is None:
        # freeze all the attributes
        for attrname, attr in self.__dict__.items():
            unparam_attr(attrname, attr)
    else:
        attr = self.get_parameter(attrname)
        if index is None:
            # freeze a scalar attribute
            unparam_attr(attrname, attr)
        else:
            # freeze an item in an array attribute
            item = attr[index]
            if isinstance(item, ParamDef):
                attr[index] = item.value


def _setattr(self, attrname: str, value: Any) -> None:
    conversion = self._conversions.get(attrname, _nop)
    try:
        # Try to convert the value
        if isinstance(value, _ACCEPTED):
            value.set_conversion(conversion)
        else:
            value = conversion(value)
    except Exception as exc:
        # Conversion failed
        exc.args = (f"{self._ident(attrname)}: {exc}",)
        raise
    else:
        # Conversion succeeded
       object.__setattr__(self, attrname, value)


def _getattribute(self, attrname: str) -> Any:
    attr = object.__getattribute__(self, attrname)
    if isinstance(attr, (ParamDef, ParamArray)):
        return attr.value
    else:
        return attr


Element.__setattr__ = _setattr
Element.__getattribute__ = _getattribute
Element.set_parameter = set_parameter
Element.is_parameterised = is_parameterised
Element.parameterise = parameterise
Element.unparameterise = unparameterise
