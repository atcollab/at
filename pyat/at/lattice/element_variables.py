"""Variables are **references** to scalar attributes of lattice elements. There are 2
kinds of element variables:

- an :py:class:`ElementVariable` is associated to an element object, and acts on all
  occurences of this object. But it will not affect any copy, neither shallow nor deep,
  of the original object,
- a :py:class:`RefptsVariable` is not associated to an element object, but to an element
  location in a :py:class:`.Lattice`. It acts on any copy of the initial lattice. A
  *ring* argument must be provided to the *set* and *get* methods to identify the
  lattice, which may be a possibly modified copy of the original lattice
"""
from __future__ import annotations
from collections.abc import Sequence
from typing import Union, Optional
import numpy as np
from .utils import Refpts, getval, setval
from .elements import Element
from .lattice_object import Lattice
from .variables import Variable

__all__ = ["RefptsVariable", "ElementVariable"]


class RefptsVariable(Variable):
    r"""A reference to a scalar attribute of :py:class:`.Lattice` elements.

    It can refer to:

    * a scalar attribute or
    * an element of an array attribute

    of one or several :py:class:`.Element`\ s of a lattice.

    A :py:class:`RefptsVariable` is not associated to element objets themselves, but
    to the location of these elements in a lattice. So a :py:class:`RefptsVariable`
    will act equally on any copy of the initial ring.
    As a consequence, a *ring* keyword argument (:py:class:`.Lattice` object) must be
    supplied for getting or setting the variable.
    """

    def __init__(
        self, refpts: Refpts, attrname: str, index: Optional[int] = None, **kwargs
    ):
        r"""
        Parameters:
            refpts:     Location of variable :py:class:`.Element`\ s
            attrname:   Attribute name
            index:      Index in the attribute array. Use :py:obj:`None` for
              scalar attributes

        Keyword Args:
            name (str):     Name of the Variable. Default: ``''``
            bounds (tuple[float, float]):   Lower and upper bounds of the
              variable value. Default: (-inf, inf)
            delta (float):  Step. Default: 1.0
        """
        self._getf = getval(attrname, index=index)
        self._setf = setval(attrname, index=index)
        self.refpts = refpts
        super().__init__(**kwargs)

    def _setfun(self, value: float, ring: Lattice = None):
        if ring is None:
            raise ValueError("Can't set values if ring is None")
        for elem in ring.select(self.refpts):
            self._setf(elem, value)

    def _getfun(self, ring: Lattice = None) -> float:
        if ring is None:
            raise ValueError("Can't get values if ring is None")
        values = np.array([self._getf(elem) for elem in ring.select(self.refpts)])
        return np.average(values)


class ElementVariable(Variable):
    r"""A reference to a scalar attribute of :py:class:`.Lattice` elements.

    It can refer to:

    * a scalar attribute or
    * an element of an array attribute

    of one or several :py:class:`.Element`\ s of a lattice.

    An :py:class:`ElementVariable` is associated to an element object, and acts on all
    occurrences of this object. But it will not affect any copy, neither shallow nor
    deep, of the original object.
    """

    def __init__(
        self,
        elements: Union[Element, Sequence[Element]],
        attrname: str,
        index: Optional[int] = None,
        **kwargs,
    ):
        r"""
        Parameters:
            elements:   :py:class:`.Element` or Sequence[Element] whose
              attribute is varied
            attrname:   Attribute name
            index:      Index in the attribute array. Use :py:obj:`None` for
              scalar attributes

        Keyword Args:
            name (str):     Name of the Variable. Default: ``''``
            bounds (tuple[float, float]):   Lower and upper bounds of the
              variable value. Default: (-inf, inf)
            delta (float):  Step. Default: 1.0
        """
        # Ensure the uniqueness of elements
        if isinstance(elements, Element):
            self.elements = {elements}
        else:
            self.elements = set(elements)
        # Check that the attribute is not already a parameter
        if any(el.is_parametrised(attrname, index=index) for el in self.elements):
            raise TypeError(f"{attrname} attribute is a Variable")
        self._getf = getval(attrname, index=index)
        self._setf = setval(attrname, index=index)
        super().__init__(**kwargs)
        self._history.append(self._getfun())

    def _setfun(self, value: float, **kwargs):
        for elem in self.elements:
            self._setf(elem, value)

    def _getfun(self, **kwargs) -> float:
        values = np.array([self._getf(elem) for elem in self.elements])
        return np.average(values)
