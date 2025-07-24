from __future__ import annotations

from typing import Any
import numpy as np
from .utils import getval, Refpts, All
from .lattice_object import Lattice
from .parameters import Param


def set_parameter(
    self, refpts: Refpts, attrname: str, value: Any, index: int | None = None
) -> None:
    """Set a parameter as an attribute of the selected elements

    Args:
        refpts:     Element selector
        attrname:   Attribute name
        value:      Parameter or value to be set
        index:      Index into an array attribute. If *value* is a
          parameter, the attribute is converted to a
          :py:class:`.ParamArray`.
    """
    for elem in self.select(refpts):
        elem.set_parameter(attrname, value, index=index)


def parameterise(
    self, refpts: Refpts, attrname: str, index: int | None = None, name: str = ""
) -> Param:
    """Convert an attribute of the selected elements into a parameter

    A single parameter is created and assigned to all the selected
    elements. Its initial value is the mean of the original values.

    Args:
        refpts:     Element selector. The default value is `All`, which includes
          all elements.
        attrname:   Name of the attribute to parameterise.
        index:      Index into an array attribute. If *value* is a
          parameter, the attribute is converted to a
          :py:class:`.ParamArray`.
        name:       Name of the created parameter

    Returns:
        param:      The created parameter
    """
    elems = self[refpts]
    getf = getval(attrname, index=index)
    vals = np.array([getf(elem) for elem in elems])
    attr = Param(np.mean(vals), name=name)
    for elem in elems:
        elem.set_parameter(attrname, attr, index=index)
    return attr


def is_parameterised(
    self, refpts: Refpts = All, attrname: str | None = None, index: int | None = None
) -> bool:
    """Checks if any of the selected elements is parameterised.

    This method evaluates whether any of the elements within the specified
    reference points are parameterised based on the provided attribute name
    and index. It iterates over the selected elements and checks their
    parameterisation status. The function returns True if at least one element
    is parameterised; otherwise, it returns False.

    Args:
        refpts:     Element selector. The default value is `All`, which includes
          all elements.
        attrname:   Name of the attribute to check. If :py:obj:`None`, all the
          attributes are checked.
        index:      Index in an array. If :py:obj:`None`, the whole
          attribute is checked

    Returns:
        bool: True if any of the selected elements is parameterised, False
        otherwise.
    """
    for elem in self.select(refpts):
        if elem.is_parameterised(attrname=attrname, index=index):
            return True
    return False


def unparameterise(
    self, refpts: Refpts = All, attrname: str | None = None, index: int | None = None
) -> None:
    """Freeze the value of attributes of the selected elements

    Args:
        refpts:     Element selector. The default value is `All`, which includes
          all elements.
        attrname:   Attribute name. If :py:obj:`None`, all the attributes
          are frozen
        index:      Index in an array. If :py:obj:`None`, the whole
          attribute is frozen
    """
    for elem in self.select(refpts):
        elem.unparameterise(attrname=attrname, index=index)


Lattice.set_parameter = set_parameter
Lattice.parameterise = parameterise
Lattice.is_parameterised = is_parameterised
Lattice.unparameterise = unparameterise
