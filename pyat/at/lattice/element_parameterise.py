from __future__ import annotations

__all__ = [
    "set_parameter",
    "get_parameter",
    "is_parameterised",
    "parameterise",
    "unparameterise",
]

from typing import Any

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
            raise IndexError(
                f"Index {index} out of bounds for {self.FamName}.{attrname}"
            ) from exc

    if index is None:
        setattr(self, attrname, value)
    else:
        array = self._get_attribute(attrname)
        if not isinstance(array, ParamArray) and isinstance(value, ParamDef):
            array = ParamArray(array, shape=array.shape, dtype=array.dtype)
            setitem(array, index, value)
            setattr(self, attrname, array)
        else:
            setitem(array, index, value)


def _get_attribute(self, attrname: str, index: int | None = None) -> Any:
    try:
        attr = self.__dict__[attrname]
    except KeyError:
        try:
            attr = self._parameters[attrname]
        except KeyError:
            raise AttributeError(
                f"{self.FamName} has no attribute '{attrname}'"
            ) from None
    if index is not None:
        try:
            attr = attr[index]
        except IndexError as exc:
            raise IndexError(f"{self.FamName}.{attrname}: {exc}") from None
    return attr


def get_parameter(self, attrname: str, index: int | None = None) -> ParamDef:
    """Extract a parameter of an element

    Unlike :py:func:`getattr`, :py:func:`get_parameter` returns the
    parameter itself instead of its value. It the item is not a parameter,
    both functions are equivalent, the value is returned. Properties cannot
    be accessed, one must use the associated array item.

    Args:
        attrname:   Attribute name
        index:      Index in an array attribute. If :py:obj:`None`, the
          whole attribute is set

    Returns:
        The parameter object

    Raises:
        TypeError if the attribute is not a Parameter
    """
    attr = self._get_attribute(attrname, index=index)
    if not isinstance(attr, ParamDef):
        idx = "" if index is None else f"[{index}]"
        message = f"\n\n{self.FamName}.{attrname}{idx} is not a parameter.\n"
        raise TypeError(message)
    return attr


def is_parameterised(
    self, attrname: str | None = None, index: int | None = None
) -> bool:
    """Check for the parametrisation of an element

    Args:
        attrname:   Attribute name. If :py:obj:`None`, return :py:obj:`True`
          if any attribute is parametrized
        index:      Index in an array attribute. If :py:obj:`None`, the
          whole attribute is tested for parametrisation

    Returns:
        True if the element, attribute, or array item is parameterised, False otherwise
    """
    if attrname is None:
        return len(self._parameters) > 0
    else:
        attr = self._get_attribute(attrname, index=index)
        return isinstance(attr, (ParamDef, ParamArray))


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
    vini = self._get_attribute(attrname, index=index)

    if isinstance(vini, ParamBase):
        return vini

    try:
        param = Param(vini, name=name)
    except TypeError as exc:
        raise TypeError(
            f"Cannot parameterise {self.FamName}.{attrname}: {exc}"
        ) from None

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
    if attrname is None:
        # freeze all the attributes
        for attrname, attr in self._parameters.copy().items():
            setattr(self, attrname, attr.value)
    else:
        attr = self._get_attribute(attrname)
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


def _setattr(self, key: str, value: Any) -> None:
    try:
        # Try to convert the value
        if isinstance(value, _ACCEPTED):
            value.set_conversion(self._conversions.get(key, _nop))
        elif not isinstance(value, ParamArray):
            value = self._conversions.get(key, _nop)(value)
    except Exception as exc:
        # Conversion failed - provide a more descriptive error message
        exc.args = (f"In element {self.FamName}, parameter {key}: {exc}",)
        raise
    else:
        # Conversion succeeded
        if isinstance(value, (ParamDef, ParamArray)):
            # Store the parameter and remove the attribute
            self._parameters[key] = value
            try:
                delattr(self, key)
            except AttributeError:
                # Attribute doesn't exist, which is fine
                pass
        else:
            # Store the attribute and remove the parameter
            object.__setattr__(self, key, value)
            try:
                del self._parameters[key]
            except KeyError:
                # Parameter doesn't exist, which is fine
                pass


def _getattr(self, key) -> Any:
    try:
        return self._parameters[key].value
    except KeyError:
        raise AttributeError(key) from None


Element.__setattr__ = _setattr
Element.__getattr__ = _getattr
Element.set_parameter = set_parameter
Element.get_parameter = get_parameter
Element.is_parameterised = is_parameterised
Element.parameterise = parameterise
Element.unparameterise = unparameterise
Element._get_attribute = _get_attribute
