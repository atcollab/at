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

__all__ = ["RefptsVariable", "ElementVariable"]

from collections.abc import Sequence, Callable
from typing import Any

import numpy as np

from .elements import Element
from .lattice_object import Lattice
from .utils import Refpts, getval, setval
from .variables import VariableBase, Number


class RefptsVariable(VariableBase):
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

    _getf: Callable[[Element], Number]
    _setf: Callable[[Element, Number], None]
    refpts: Refpts

    def __init__(
        self, refpts: Refpts, attrname: str, index: int | None = None, **kwargs
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
            ring (Lattice): Default lattice. It will be used if *ring* is not provided
              to the :py:meth:`~.VariableBase.get` or :py:meth:`~.VariableBase.set`
              methods. Note that if *ring* is fixed, you should consider using a
              :py:class:`.ElementVariable` instead.
        """
        self._getf = getval(attrname, index=index)
        self._setf = setval(attrname, index=index)
        self.refpts = refpts
        super().__init__(**kwargs)

    def _setfun(self, value: Number, ring: Lattice | None = None, **_) -> None:
        if ring is None:
            raise ValueError("Can't set values if ring is None.\n"
                             "Try to use an ElementVariable if possible")
        for elem in ring.select(self.refpts):
            self._setf(elem, value)

    def _getfun(self, ring: Lattice | None = None, **_) -> Number:
        if ring is None:
            raise ValueError("Can't get values if ring is None")
        values = np.array([self._getf(elem) for elem in ring.select(self.refpts)])
        return np.average(values)


class ElementVariable(VariableBase):
    r"""A reference to a scalar attribute of :py:class:`.Lattice` elements.

    It can refer to:

    * a scalar attribute or
    * an element of an array attribute

    of one or several :py:class:`.Element`\ s of a lattice.

    An :py:class:`ElementVariable` is associated to an element object, and acts on all
    occurrences of this object. But it will not affect any copy, neither shallow nor
    deep, of the original object.
    """

    _getf: Callable[[Element], Any]
    _setf: Callable[[Element, Any], None]
    _elements: set[Element]

    def __init__(
        self,
        elements: Element | Sequence[Element],
        attrname: str,
        index: int | None = None,
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
            self._elements = {elements}
        else:
            self._elements = set(elements)
        self._getf = getval(attrname, index=index)
        self._setf = setval(attrname, index=index)
        super().__init__(**kwargs)

    def _setfun(self, value: Number, **_) -> None:
        for elem in self._elements:
            self._setf(elem, value)

    def _getfun(self, **_) -> Number:
        values = np.array([self._getf(elem) for elem in self._elements])
        return np.average(values)

    @property
    def elements(self) -> set[Element]:
        """Return the set of elements acted upon by the variable"""
        return self._elements
