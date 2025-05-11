from __future__ import annotations

__all__ = ["set_parameter", "get_parameter", "is_parametrised", "parametrise"]

import numpy as np

from .elements import Element
from .variables import _nop
from .parameters import Param, ParamBase, ParamArray


def set_parameter(self, attrname: str, value, index: int | None = None) -> None:
    """Set an element's parameter

    This allows setting a parameter into an item of an array attribute.

    Args:
        attrname:   Attribute name
        value:      Parameter or value to be set
        index:      Index into an array attribute. If *value* is a
          parameter, the attribute is converted to a
          :py:class:`.ParamArray`.
    """
    if index is None:
        setattr(self, attrname, value)
    else:
        attr = self._get_attribute(attrname)
        attr[index] = value


def _get_attribute(self, attrname: str, index: int | None = None):
    attr = self.__dict__[attrname]
    if index is not None:
        attr = attr[index]
    return attr


def get_parameter(self, attrname: str, index: int | None = None):
    """Extract a parameter of an element

    Unlike :py:func:`getattr`, :py:func:`get_parameter` returns the
    parameter itself instead of its value. It the item is not a parameter,
    both functions are equivalent, the value is returned. Properties cannot
    be accessed, one must use the associated array item.

    Args:
        attrname:   Attribute name
        index:      Index in an array attribute. If :py:obj:`None`, the
          whole attribute is set
    """
    attr = self._get_attribute(attrname, index=index)
    if not isinstance(attr, ParamBase):
        idx = "" if index is None else f"[{index}]"
        message = f"\n\n{self.FamName}.{attrname}{idx} is not a parameter.\n"
        raise TypeError(message)
    return attr


def is_parametrised(self, attrname: str | None = None,
                    index: int | None = None) -> bool:
    """Check for the parametrisation of an element

    Args:
        attrname:   Attribute name. If :py:obj:`None`, return :py:obj:`True`
          if any attribute is parametrized
        index:      Index in an array attribute. If :py:obj:`None`, the
          whole attribute is tested for parametrisation
    """
    if attrname is None:
        for attr in self.__dict__:
            if self.is_parametrised(attr):
                return True
        return False
    else:
        attr = self._get_attribute(attrname, index=index)
        if isinstance(attr, ParamBase):
            return True
        elif isinstance(attr, np.ndarray):
            return any(isinstance(item, ParamBase) for item in attr.flat)
        else:
            return False


def parametrise(self, attrname: str, index: int | None = None,
                name: str = '') -> ParamBase:
    """Convert an attribute into a parameter

    The value of the attribute is kept unchanged. If the attribute is
    already parametrised, the existing parameter is returned.

    Args:
        attrname:   Attribute name
        index:      Index in an array. If :py:obj:`None`, the
          whole attribute is parametrised
        name:       Name of the created parameter

    Returns:
        param:      A :py:class:`.ParamArray` for an array attribute,
          a :py:class:`.Param` for a scalar attribute or an item in an
          array attribute

    """
    vini = self._get_attribute(attrname, index=index)

    if isinstance(vini, ParamBase):
        return vini

    attr = Param(vini, name=name)   # raises TypeError if vini is not a Number

    if index is None:
        setattr(self, attrname, attr)
    else:
        varr = self._get_attribute(attrname)
        varr[index] = attr          # raises IndexError it the attr is not an array
    return attr


def unparametrise(self, attrname: str | None = None,
                  index: int | None = None) -> None:
    """Freeze the parameter values

    Args:
        attrname:   Attribute name. If :py:obj:`None`, all the attributes
          are frozen
        index:      Index in an array. If :py:obj:`None`, the whole
          attribute is frozen

    Attributes which are not parameters are silently ignored.
    """

    def unparam_attr(attrname, attr):
        if isinstance(attr, ParamBase):
            setattr(self, attrname, attr.value)
        elif isinstance(attr, np.ndarray):
            for i, item in enumerate(attr.flat):
                if isinstance(item, ParamBase):
                    ij = np.unravel_index(i, attr.shape)
                    attr[ij] = item.value

    if attrname is None:
        for attrname, attr in self.__dict__.items():
            unparam_attr(attrname, attr)
    else:
        attr = self._get_attribute(attrname)
        if index is None:
            unparam_attr(attrname, attr)
        else:
            item = attr[index]
            if isinstance(item, ParamBase):
                attr[index] = item.value


def _setattr(self, key, value):
    try:
        if isinstance(value, ParamBase):
            value.set_conversion(self._conversions.get(key, _nop))
        else:
            value = self._conversions.get(key, _nop)(value)
    except Exception as exc:
        exc.args = ('In element {0}, parameter {1}: {2}'.format(
            self.FamName, key, exc),)
        raise
    else:
        super(Element, self).__setattr__(key, value)


def _getattribute(self, key):
    attr = super(Element, self).__getattribute__(key)
    if isinstance(attr, (ParamBase, ParamArray)):
        return attr.value
    else:
        return attr

Element.__setattr__ = _setattr
Element.__getattribute__ = _getattribute
Element.set_parameter = set_parameter
Element.get_parameter = get_parameter
Element.is_parametrised = is_parametrised
Element.parametrise = parametrise
Element.unparametrise = unparametrise
Element._get_attribute = _get_attribute
