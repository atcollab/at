r"""
Helper functions for working with AT lattices.

A lattice as understood by pyAT is a sequence of elements.  These functions
are useful for working with these sequences.

.. _refpts:

**Selecting elements in a lattice:**

The *refpts* argument allows functions to select locations in the lattice. The
location is defined as the entrance of the selected element. *refpts* may be:

#. an integer in the range *[-len(lattice), len(lattice)-1]*
   selecting an element according to the python indexing rules.
   As a special case, *len(lattice)* is allowed and refers to the end
   of the last element,
#. an ordered list of such integers without duplicates,
#. a numpy array of booleans of maximum length *len(lattice)+1*,
   where selected elements are :py:obj:`True`,
#. :py:obj:`None`, meaning an empty selection,
#. :py:obj:`.All`, meaning "all possible reference points": the entrance of all
   elements plus the end of the last element,
#. :py:obj:`.End`, selecting the end of the last element,
#. an element type, selecting all the elements of that type in
   the lattice, e.g. :pycode:`at.Sextupole`,
#. a string, selecting all the elements whose `FamName` attribute matches it.
   Unix shell-style wildcards are accepted, e.g. `"Q[FD]*"`,
#. a callable :pycode:`filtfunc` such that :pycode:`filtfunc(elem)`
   is :py:obj:`True` for selected elements.

"""

from __future__ import annotations

import functools
import re
from enum import Enum
from fnmatch import fnmatch
from itertools import compress
from operator import attrgetter
from typing import Union

# Necessary for type aliases in python <= 3.8 :
# from collections.abc import Callable, Sequence, Iterator
from typing import Callable, Sequence, Iterator, Type

import numpy
import numpy.typing as npt

from .elements import Element, Dipole

_GEOMETRY_EPSIL = 1.0e-3

__all__ = [
    "All",
    "End",
    "AtError",
    "AtWarning",
    "BoolRefpts",
    "Uint32Refpts",
    "check_radiation",
    "check_6d",
    "set_radiation",
    "set_6d",
    "make_copy",
    "uint32_refpts",
    "bool_refpts",
    "get_uint32_index",
    "get_bool_index",
    "checkattr",
    "checktype",
    "checkname",
    "get_elements",
    "get_s_pos",
    "refpts_count",
    "refpts_iterator",
    "get_value_refpts",
    "set_value_refpts",
    "Refpts",
    "get_geometry",
    "setval",
    "getval",
]


class AtError(Exception):
    pass


class AtWarning(UserWarning):
    pass


_typ1 = "None, All, End, int, bool"

_typ2 = "None, All, End, int, bool, str, type[Element], ElementFilter"


class RefptsCode(Enum):
    All = "All"
    End = "End"


ElementFilter = Callable[[Element], bool]
BoolRefpts = npt.NDArray[bool]
Uint32Refpts = npt.NDArray[numpy.uint32]
RefIndex = Union[None, int, Sequence[int], bool, Sequence[bool], RefptsCode]
Refpts = Union[Type[Element], Element, ElementFilter, str, RefIndex]


#: :py:obj:`All` is a special value to be used as *refpts*. It means
#: "all possible reference points": the entrance of all elements plus the end
#: of the last element.
All = RefptsCode.All

#: :py:obj:`End` is a special value to be used as *refpts*. It refers to the
#: end of the last element.
End = RefptsCode.End


def _chkattr(attrname: str, el):
    return hasattr(el, attrname)


def _chkattrval(attrname: str, attrvalue, el):
    try:
        v = getattr(el, attrname)
    except AttributeError:
        return False
    else:
        return v == attrvalue


def _chkpattern(pattern: str, el):
    return fnmatch(el.FamName, pattern)


def _chkregex(pattern: str, el):
    rgx = re.compile(pattern)
    return rgx.fullmatch(el.FamName)


def _chktype(eltype: type, el):
    return isinstance(el, eltype)


def _type_error(refpts, types):
    if isinstance(refpts, numpy.ndarray):
        tp = refpts.dtype.type
    else:
        tp = type(refpts)
    return TypeError(f"Invalid refpts type {tp}. Allowed types: {types}")


# setval and getval return pickleable functions: no inner, nested function
# is allowed. So nested functions are replaced be module-level callable
# class instances
class _AttrItemGetter:
    __slots__ = ["attrname", "index"]

    def __init__(self, attrname: str, index: int):
        self.attrname = attrname
        self.index = index

    def __call__(self, elem):
        return getattr(elem, self.attrname)[self.index]


def getval(attrname: str, index: int | None = None) -> Callable:
    """Return a callable object which fetches item *index* of
    attribute *attrname* of its operand. Examples:

    - After ``f = getval('Length')``, ``f(elem)`` returns ``elem.Length``
    - After ``f = getval('PolynomB', index=1)``, ``f(elem)`` returns
      ``elem.PolynomB[1]``

    """
    if index is None:
        return attrgetter(attrname)
    else:
        return _AttrItemGetter(attrname, index)


class _AttrSetter:
    __slots__ = ["attrname"]

    def __init__(self, attrname: str):
        self.attrname = attrname

    def __call__(self, elem, value):
        setattr(elem, self.attrname, value)


class _AttrItemSetter:
    __slots__ = ["attrname", "index"]

    def __init__(self, attrname: str, index: int):
        self.attrname = attrname
        self.index = index

    def __call__(self, elem, value):
        getattr(elem, self.attrname)[self.index] = value


def setval(attrname: str, index: int | None = None) -> Callable:
    """Return a callable object which sets the value of  item *index* of
    attribute *attrname* of its 1st argument to it 2nd orgument.

    - After ``f = setval('Length')``, ``f(elem, value)`` is equivalent to
      ``elem.Length = value``
    - After ``f = setval('PolynomB', index=1)``, ``f(elem, value)`` is
      equivalent to ``elem.PolynomB[1] = value``

    """
    if index is None:
        return _AttrSetter(attrname)
    else:
        return _AttrItemSetter(attrname, index)


def check_radiation(rad: bool) -> Callable:
    r"""Deprecated decorator for optics functions (see :py:func:`check_6d`).

    :meta private:
    """
    return check_6d(rad)


def set_radiation(rad: bool) -> Callable:
    r"""Deprecated decorator for optics functions (see :py:func:`set_6d`).

    :meta private:
    """
    return set_6d(rad)


def check_6d(is_6d: bool) -> Callable:
    r"""Decorator for optics functions

    Wraps a function like :pycode:`func(ring, *args, **kwargs)` where
    :pycode:`ring` is a :py:class:`.Lattice` object. Raise :py:class:`.AtError`
    if :pycode:`ring.is_6d` is not the desired *is_6d* state. No test is
    performed if :pycode:`ring` is not a :py:class:`.Lattice`.

    Parameters:
        is_6d: Desired 6D state

    Raises:
        AtError: if :pycode:`ring.is_6d` is not *is_6d*

    Example:

        >>> @check_6d(False)
        ... def find_orbit4(ring, dct=None, guess=None, **kwargs):
                ...

        Raises :py:class:`.AtError` if :pycode:`ring.is_6d` is :py:obj:`True`

    See Also:
        :py:func:`set_6d`
    """

    def radiation_decorator(func):
        @functools.wraps(func)
        def wrapper(ring, *args, **kwargs):
            ringrad = getattr(ring, "is_6d", is_6d)
            if ringrad != is_6d:
                raise AtError(f'{func.__name__} needs "ring.is_6d" {is_6d}')
            return func(ring, *args, **kwargs)

        return wrapper

    return radiation_decorator


def set_6d(is_6d: bool) -> Callable:
    # noinspection PyUnresolvedReferences
    r"""Decorator for optics functions

    Wraps a function like :pycode:`func(ring, *args, **kwargs)` where
    :pycode:`ring` is a :py:class:`.Lattice` object. If :pycode:`ring.is_6d`
    is not the desired *is_6d* state, :pycode:`func()` will be called with a
    modified copy of :pycode:`ring` satisfying the *is_6d* state.

    Parameters:
        is_6d: Desired 6D state

    Example:

        >>> @set_6d(True)
        ... def compute(ring)
                ...

        is roughly equivalent to:

        >>> if ring.is_6d:
        ...     compute(ring)
        ... else:
        ...     compute(ring.enable_6d(copy=True))

        See Also:
            :py:func:`check_6d`, :py:meth:`.Lattice.enable_6d`,
            :py:meth:`.Lattice.disable_6d`
    """
    if is_6d:

        def setrad_decorator(func):
            @functools.wraps(func)
            def wrapper(ring, *args, **kwargs):
                rg = ring if ring.is_6d else ring.enable_6d(copy=True)
                return func(rg, *args, **kwargs)

            return wrapper

    else:

        def setrad_decorator(func):
            @functools.wraps(func)
            def wrapper(ring, *args, **kwargs):
                rg = ring.disable_6d(copy=True) if ring.is_6d else ring
                return func(rg, *args, **kwargs)

            return wrapper

    return setrad_decorator


def make_copy(copy: bool) -> Callable:
    r"""Decorator for optics functions

    Wraps a function like :pycode:`func(ring, refpts, *args, **kwargs)` where
    :pycode:`ring` is a :py:class:`.Lattice` object.

    If *copy* is :py:obj:`False`, the function is not modified,

    If *copy* is :py:obj:`True`, a shallow copy of :pycode:`ring` is done,
    then the elements selected by *refpts* are deep-copied, then :pycode:`func`
    is applied to the copy, and the new :pycode:`ring` is returned.

    Parameters:
        copy:   If :py:obj:`True`, apply the decorated function to a copy of
          :pycode:`ring`
    """
    if copy:

        def copy_decorator(func):
            @functools.wraps(func)
            def wrapper(ring, refpts, *args, **kwargs):
                try:
                    ring = ring.replace(refpts)
                except AttributeError:
                    check = get_bool_index(ring, refpts)
                    ring = [el.deepcopy() if ok else el for el, ok in zip(ring, check)]
                func(ring, refpts, *args, **kwargs)
                return ring

            return wrapper

    else:

        def copy_decorator(func):
            return func

    return copy_decorator


def uint32_refpts(
    refpts: RefIndex, n_elements: int, endpoint: bool = True, types: str = _typ1
) -> Uint32Refpts:
    r"""Return a :py:obj:`~numpy.uint32` array of element indices selecting
    ring elements.

    Parameters:
        refpts:     Element selector. *refpts* may be:

          #. an integer or a sequence of integers
             (0 indicating the first element),
          #. a sequence of booleans marking the selected elements,
          #. :py:obj:`None`, meaning empty selection,
          #. :py:obj:`.All`, meaning "all possible reference points",
          #. :py:obj:`.End`, selecting the end of the last element.
        n_elements: Length of the sequence of elements
        endpoint:   if :py:obj:`True`, allow *n_elements* as a
          special index, referring to the end of the last element.
        types:      Allowed types

    Returns:
        uint32_ref (Uint32Refpts):  :py:obj:`~numpy.uint32` numpy array used
          for indexing :py:class:`.Element`\ s in a lattice.
    """
    refs = numpy.ravel(refpts)
    if refpts is RefptsCode.All:
        stop = n_elements + 1 if endpoint else n_elements
        return numpy.arange(stop, dtype=numpy.uint32)
    elif refpts is RefptsCode.End:
        if not endpoint:
            raise IndexError('"End" index out of range')
        return numpy.array([n_elements], dtype=numpy.uint32)
    elif (refpts is None) or (refs.size == 0):
        return numpy.array([], dtype=numpy.uint32)
    elif numpy.issubdtype(refs.dtype, numpy.bool_):
        return numpy.flatnonzero(refs).astype(numpy.uint32)
    elif numpy.issubdtype(refs.dtype, numpy.integer):
        # Handle negative indices
        if endpoint:
            refs = numpy.array(
                [i if (i == n_elements) else i % n_elements for i in refs],
                dtype=numpy.uint32,
            )
        else:
            refs = numpy.array([i % n_elements for i in refs], dtype=numpy.uint32)
        # Check ascending
        if refs.size > 1:
            prev = refs[0]
            for nxt in refs[1:]:
                if nxt < prev:
                    raise IndexError("Index out of range or not in ascending order")
                elif nxt == prev:
                    raise IndexError("Duplicated index")
                prev = nxt

        return refs
    else:
        raise _type_error(refpts, types)


# noinspection PyIncorrectDocstring
def get_uint32_index(
    ring: Sequence[Element], refpts: Refpts, endpoint: bool = True, regex: bool = False
) -> Uint32Refpts:
    # noinspection PyUnresolvedReferences, PyShadowingNames
    r"""Returns an integer array of element indices, selecting ring elements.

    Parameters:
        refpts:         Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        endpoint:   if :py:obj:`True`, allow *len(ring)* as a
          special index, referring to the end of the last element.
        regex: Use regular expression for *refpts* string matching instead of
          Unix shell-style wildcards.

    Returns:
        uint32_ref (Uint32Refpts): uint32 numpy array used for indexing
          :py:class:`.Element`\ s in a lattice.

    Examples:

        >>> get_uint32_index(ring, at.Sextupole)
        array([ 21,  27,  35,  89,  97, 103], dtype=uint32)

        numpy array of indices of all :py:class:`.Sextupole`\ s

        >>> get_uint32_index(ring, at.End)
        array([121], dtype=uint32)

        numpy array([:pycode:`len(ring)+1`])

        >>> get_uint32_index(ring, at.checkattr("Frequency"))
        array([0], dtype=uint32)

        numpy array of indices of all elements having a 'Frequency'
        attribute
    """
    if isinstance(refpts, type):
        checkfun = checktype(refpts)
    elif callable(refpts):
        checkfun = refpts
    elif isinstance(refpts, Element):
        checkfun = checktype(type(refpts))
    elif isinstance(refpts, str):
        checkfun = checkname(refpts, regex=regex)
    else:
        return uint32_refpts(refpts, len(ring), endpoint=endpoint, types=_typ2)

    return numpy.fromiter(
        (i for i, el in enumerate(ring) if checkfun(el)), dtype=numpy.uint32
    )


def bool_refpts(
    refpts: RefIndex, n_elements: int, endpoint: bool = True, types: str = _typ1
) -> BoolRefpts:
    r"""Returns a :py:class:`bool` array of element indices, selecting ring
    elements.

    Parameters:
        refpts:     Element selector. *refpts* may be:

          #. an integer or a sequence of integers
             (0 indicating the first element),
          #. a sequence of booleans marking the selected elements,
          #. :py:obj:`None`, meaning empty selection,
          #. :py:obj:`.All`, meaning "all possible reference points",
          #. :py:obj:`.End`, selecting the end of the last element.
        n_elements: Length of the lattice
        endpoint:   if :py:obj:`True`, allow *n_elements* as a
          special index, referring to the end of the last element.
        types:      Allowed types

    Returns:
        bool_refs (BoolRefpts):  A bool numpy array used for indexing
          :py:class:`.Element`\ s in a lattice.
    """
    refs = numpy.ravel(refpts)
    stop = n_elements + 1 if endpoint else n_elements
    if refpts is RefptsCode.All:
        return numpy.ones(stop, dtype=bool)
    elif refpts is RefptsCode.End:
        if not endpoint:
            raise IndexError('"End" index out of range')
        brefpts = numpy.zeros(stop, dtype=bool)
        brefpts[n_elements] = True
        return brefpts
    elif (refpts is None) or (refs.size == 0):
        return numpy.zeros(stop, dtype=bool)
    elif numpy.issubdtype(refs.dtype, numpy.bool_):
        diff = stop - refs.size
        if diff <= 0:
            return refs[:stop]
        else:
            return numpy.append(refs, numpy.zeros(diff, dtype=bool))
    elif numpy.issubdtype(refs.dtype, numpy.integer):
        brefpts = numpy.zeros(stop, dtype=bool)
        brefpts[refs] = True
        return brefpts
    else:
        raise _type_error(refpts, types)


# noinspection PyIncorrectDocstring
def get_bool_index(
    ring: Sequence[Element], refpts: Refpts, endpoint: bool = True, regex: bool = False
) -> BoolRefpts:
    # noinspection PyUnresolvedReferences, PyShadowingNames
    r"""Returns a bool array of element indices, selecting ring elements.

    Parameters:
        refpts:         Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        endpoint:   if :py:obj:`True`, allow *len(ring)* as a
          special index, referring to the end of the last element.
        regex: Use regular expression for *refpts* string matching instead of
          Unix shell-style wildcards.

    Returns:
        bool_refs (BoolRefpts):  A bool numpy array used for indexing
          :py:class:`.Element`\ s in a lattice.

    Examples:

        >>> refpts = get_bool_index(ring, at.Quadrupole)

        Returns a numpy array of booleans where all :py:class:`.Quadrupole`
        are :py:obj:`True`

        >>> refpts = get_bool_index(ring, "Q[FD]*")

        Returns a numpy array of booleans where all elements whose *FamName*
        matches "Q[FD]*" are :py:obj:`True`

        >>> refpts = get_bool_index(ring, at.checkattr("K", 0.0))

        Returns a numpy array of booleans where all elements whose *K*
        attribute is 0.0 are :py:obj:`True`

        >>> refpts = get_bool_index(ring, None)

        Returns a numpy array of *len(ring)+1* :py:obj:`False` values
    """
    if isinstance(refpts, type):
        checkfun = checktype(refpts)
    elif callable(refpts):
        checkfun = refpts
    elif isinstance(refpts, Element):
        checkfun = checktype(type(refpts))
    elif isinstance(refpts, str):
        checkfun = checkname(refpts, regex=regex)
    else:
        return bool_refpts(refpts, len(ring), endpoint=endpoint, types=_typ2)

    boolrefs = numpy.fromiter(map(checkfun, ring), dtype=bool, count=len(ring))
    if endpoint:
        boolrefs = numpy.append(boolrefs, False)
    return boolrefs


def checkattr(attrname: str, attrvalue=None) -> ElementFilter:
    # noinspection PyUnresolvedReferences
    r"""Checks the presence or the value of an attribute

    Returns a function to be used as an :py:class:`.Element` filter, which
    checks the presence or the value of an attribute of the
    provided :py:class:`.Element`.
    This function can be used to extract from a ring all elements
    having a given attribute.

    Parameters:
        attrname: Attribute name
        attrvalue: Attribute value. If absent, the returned function checks
          the presence of an *attrname* attribute. If present, the
          returned function checks if :pycode:`attrname == attrvalue`.

    Returns:
        checkfun (ElementFilter):   Element filter function

    Examples:

        >>> cavs = filter(checkattr("Frequency"), ring)

        Returns an iterator over all elements in *ring* that have a
        :pycode:`Frequency` attribute

        >>> elts = filter(checkattr("K", 0.0), ring)

        Returns an iterator over all elements of ring that have a
        :pycode:`K` attribute equal to 0.0
    """
    if attrvalue is None:
        return functools.partial(_chkattr, attrname)
    else:
        return functools.partial(_chkattrval, attrname, attrvalue)


def checktype(eltype: type | tuple[type, ...]) -> ElementFilter:
    # noinspection PyUnresolvedReferences
    r"""Checks the type of an element

    Returns a function to be used as an :py:class:`.Element` filter, which
    checks the type of the provided :py:class:`.Element`.
    This function can be used to extract from a ring all elements
    having a given type.

    Parameters:
        eltype: Desired :py:class:`.Element` type

    Returns:
        checkfun (ElementFilter):   Element filter function

    Examples:

        >>> qps = filter(checktype(at.Quadrupole), ring)

        Returns an iterator over all quadrupoles in ring
    """
    return functools.partial(_chktype, eltype)


def checkname(pattern: str, regex: bool = False) -> ElementFilter:
    # noinspection PyUnresolvedReferences
    r"""Checks the name of an element

    Returns a function to be used as an :py:class:`.Element` filter,
    which checks the name of the provided :py:class:`.Element`.
    This function can be used to extract from a ring all elements
    having a given name.

    Parameters:
        pattern: Desired :py:class:`.Element` name. Unix shell-style
          wildcards are supported (see :py:func:`fnmatch.fnmatch`)
        regex: Use regular expression for *refpts* string matching instead of
          Unix shell-style wildcards.

    Returns:
        checkfun (ElementFilter):   Element filter function

    Examples:

        >>> qps = filter(checkname("QF*"), ring)

        Returns an iterator over all the elements with name starting with ``QF``.
    """
    if regex:
        return functools.partial(_chkregex, pattern)
    else:
        return functools.partial(_chkpattern, pattern)


def refpts_iterator(
    ring: Sequence[Element], refpts: Refpts, regex: bool = False
) -> Iterator[Element]:
    r"""Return an iterator over selected elements in a lattice

    Parameters:
        ring:           Lattice description
        refpts:         Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        regex: Use regular expression for *refpts* string matching instead of
          Unix shell-style wildcards.

    Returns:
        elem_iter (Iterator[Element]):  Iterator over the elements in *ring*
          selected by *refpts*.
    """
    if isinstance(refpts, type):
        checkfun = checktype(refpts)
    elif callable(refpts):
        checkfun = refpts
    elif isinstance(refpts, Element):
        checkfun = checktype(type(refpts))
    elif isinstance(refpts, str):
        checkfun = checkname(refpts, regex=regex)
    else:
        refs = numpy.ravel(refpts)
        if refpts is RefptsCode.All:
            return (el for el in ring)
        elif refpts is RefptsCode.End:
            raise IndexError('"End" is not allowed for endpoint=False')
        elif (refpts is None) or (refs.size == 0):
            return iter(())
        elif numpy.issubdtype(refs.dtype, numpy.bool_):
            return compress(ring, refs)
        elif numpy.issubdtype(refs.dtype, numpy.integer):
            return (ring[i] for i in refs)
        else:
            raise _type_error(refpts, _typ2)

    return filter(checkfun, ring)


# noinspection PyUnusedLocal,PyIncorrectDocstring
def refpts_count(
    refpts: RefIndex, n_elements: int, endpoint: bool = True, types: str = _typ1
) -> int:
    r"""Returns the number of reference points

    Parameters:
        refpts:     refpts may be:

          #. an integer or a sequence of integers
             (0 indicating the first element),
          #. a sequence of booleans marking the selected elements,
          #. :py:obj:`None`, meaning empty selection,
          #. :py:obj:`.All`, meaning "all possible reference points",
          #. :py:obj:`.End`, selecting the end of the last element.
        n_elements: Lattice length
        endpoint:   if :py:obj:`True`, allow *n_elements* as a
          special index, referring to the end of the last element.

    Returns:
        nrefs (int):  The number of reference points
    """
    refs = numpy.ravel(refpts)
    if refpts is RefptsCode.All:
        return n_elements + 1 if endpoint else n_elements
    elif refpts is RefptsCode.End:
        if not endpoint:
            raise IndexError('"End" index out of range')
        return 1
    elif (refpts is None) or (refs.size == 0):
        return 0
    elif numpy.issubdtype(refs.dtype, numpy.bool_):
        return numpy.count_nonzero(refs)
    elif numpy.issubdtype(refs.dtype, numpy.integer):
        return len(refs)
    else:
        raise _type_error(refpts, types)


def _refcount(
    ring: Sequence[Element], refpts: Refpts, endpoint: bool = True, regex: bool = False
) -> int:
    # noinspection PyUnresolvedReferences, PyShadowingNames
    r"""Returns the number of reference points

    Parameters:
        refpts:         Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        endpoint:   if :py:obj:`True`, allow *len(ring)* as a
          special index, referring to the end of the last element.
        regex: Use regular expression for *refpts* string matching instead of
          Unix shell-style wildcards.

    Returns:
        nrefs (int):  The number of reference points

    Examples:

        >>> refpts = ring.refcount(at.Sextupole)
        6

        Returns the number of :py:class:`.Sextupole`\ s in the lattice

        >>> refpts = ring.refcount(at.All)
        122

        Returns *len(ring)+1*

        >>> refpts = ring.refcount(at.All, endpoint=False)
        121

        Returns *len(ring)*
    """
    if isinstance(refpts, type):
        checkfun = checktype(refpts)
    elif callable(refpts):
        checkfun = refpts
    elif isinstance(refpts, Element):
        checkfun = checktype(type(refpts))
    elif isinstance(refpts, str):
        checkfun = checkname(refpts, regex=regex)
    else:
        return refpts_count(refpts, len(ring), endpoint=endpoint, types=_typ2)

    return len(list(filter(checkfun, ring)))


# noinspection PyUnusedLocal,PyIncorrectDocstring
def get_elements(ring: Sequence[Element], refpts: Refpts, regex: bool = False) -> list:
    r"""Returns a list of elements selected by *key*.

    Deprecated: :pycode:`get_elements(ring, refpts)` is :pycode:`ring[refpts]`

    Parameters:
        ring:           Lattice description
        refpts:         Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        regex: Use regular expression for *refpts* string matching instead of
          Unix shell-style wildcards.

    Returns:
        elem_list (list):  list of :py:class:`.Element`\ s matching key
    """
    return list(refpts_iterator(ring, refpts, regex=regex))


def get_value_refpts(
    ring: Sequence[Element],
    refpts: Refpts,
    attrname: str,
    index: int | None = None,
    regex: bool = False,
):
    r"""Extracts attribute values from selected
        lattice :py:class:`.Element`\ s.

    Parameters:
        ring:           Lattice description
        refpts:         Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        attrname:   Attribute name
        index:      index of the value to retrieve if *attrname* is
          an array. If :py:obj:`None` the full array is retrieved
        regex: Use regular expression for *refpts* string matching instead of
          Unix shell-style wildcards.

    Returns:
        attrvalues: numpy Array of attribute values.
    """
    getf = getval(attrname, index=index)
    return numpy.array(
        [getf(elem) for elem in refpts_iterator(ring, refpts, regex=regex)]
    )


def set_value_refpts(
    ring: Sequence[Element],
    refpts: Refpts,
    attrname: str,
    attrvalues,
    index: int | None = None,
    increment: bool = False,
    copy: bool = False,
    regex: bool = False,
):
    r"""Set the values of an attribute of an array of elements based on
    their refpts


    Parameters:
        ring:       Lattice description
        refpts:         Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        attrname:   Attribute name
        attrvalues: Attribute values
        index:      index of the value to set if *attrname* is
          an array. if :py:obj:`None`, the full array is replaced by
          *attrvalue*
        increment:  If :py:obj:`True`, *attrvalues* are added to the
          initial values.
        regex: Use regular expression for *refpts* string matching
          instead of Unix shell-style wildcards.
        copy:       If :py:obj:`False`, the modification is done in-place,

          If :py:obj:`True`, return a shallow copy of the lattice. Only the
          modified elements are copied.

          .. Caution::

             a shallow copy means that all non-modified
             elements are shared with the original lattice.
             Any further modification will affect both lattices.
    """
    setf = setval(attrname, index=index)
    if increment:
        attrvalues += get_value_refpts(ring, refpts, attrname, index=index, regex=regex)
    else:
        attrvalues = numpy.broadcast_to(
            attrvalues, (_refcount(ring, refpts, regex=regex),)
        )

    # noinspection PyShadowingNames
    @make_copy(copy)
    def apply(ring, refpts, values, regex):
        for elm, val in zip(refpts_iterator(ring, refpts, regex=regex), values):
            setf(elm, val)

    return apply(ring, refpts, attrvalues, regex)


def get_s_pos(
    ring: Sequence[Element], refpts: Refpts = All, regex: bool = False
) -> Sequence[float]:
    # noinspection PyUnresolvedReferences
    r"""Returns the locations of selected elements

    Parameters:
        ring:       Lattice description
        refpts:     Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        regex: Use regular expression for *refpts* string matching instead of
          Unix shell-style wildcards.

    Returns:
        s_pos:  Array of locations of the elements selected by *refpts*

    Example:

        >>> get_s_pos(ring, at.End)
        array([26.37428795])

        Position at the end of the last element: length of the lattice
    """
    # Positions at the end of each element.
    s_pos = numpy.cumsum([getattr(el, "Length", 0.0) for el in ring])
    # Prepend position at the start of the first element.
    s_pos = numpy.concatenate(([0.0], s_pos))
    return s_pos[get_bool_index(ring, refpts, regex=regex)]


def get_geometry(
    ring: list[Element],
    refpts: Refpts = All,
    start_coordinates: tuple[float, float, float] = (0, 0, 0),
    centered: bool = False,
    regex: bool = False,
):
    # noinspection PyShadowingNames
    r"""Compute the 2D ring geometry in cartesian coordinates

    Parameters:
        ring:               Lattice description.
        refpts:     Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        start_coordinates:  *x*, *y*, *angle* at starting point. *x*
          and *y* are ignored if *centered* is :py:obj:`True`.
        centered:           if :py:obj:`True` the coordinates origin is the
          centre of the ring.
        regex: Use regular expression for *refpts* string matching instead of
          Unix shell-style wildcards.

    Returns:
        geomdata:           recarray containing, x, y, angle.
        radius:             machine radius at the beginning of the lattice.

            .. attention::
               This radius is different from the radius usually defined as
               :math:`C/2\pi`

    Example:

       >>> geomdata, radius = get_geometry(ring)
    """

    geom_dtype = [("x", numpy.float64), ("y", numpy.float64), ("angle", numpy.float64)]
    boolrefs = get_bool_index(ring, refpts, endpoint=True, regex=regex)
    nrefs = refpts_count(boolrefs, len(ring))
    geomdata = numpy.recarray((nrefs,), dtype=geom_dtype)
    xx = numpy.zeros(len(ring) + 1)
    yy = numpy.zeros(len(ring) + 1)
    angle = numpy.zeros(len(ring) + 1)
    x0, y0, t0 = start_coordinates
    x, y = 0.0, 0.0
    t = t0

    xx[0] = x
    yy[0] = y
    angle[0] = t
    for ind, el in enumerate(ring):
        ll = el.Length
        if isinstance(el, Dipole) and el.BendingAngle != 0:
            ang = 0.5 * el.BendingAngle
            ll *= numpy.sin(ang) / ang
        else:
            ang = 0.0
        t -= ang
        x += ll * numpy.cos(t)
        y += ll * numpy.sin(t)
        t -= ang
        xx[ind + 1] = x
        yy[ind + 1] = y
        angle[ind + 1] = t

    dff = (t + _GEOMETRY_EPSIL) % (2.0 * numpy.pi) - _GEOMETRY_EPSIL
    if abs(dff) < _GEOMETRY_EPSIL:
        xcenter = numpy.mean(xx)
        ycenter = numpy.mean(yy)
    elif abs(dff - numpy.pi) < _GEOMETRY_EPSIL:
        xcenter = 0.5 * x
        ycenter = 0.5 * y
    else:
        num = numpy.cos(t) * x + numpy.sin(t) * y
        den = numpy.sin(t - t0)
        xcenter = -num * numpy.sin(t0) / den
        ycenter = num * numpy.cos(t0) / den
    radius = numpy.sqrt(xcenter * xcenter + ycenter * ycenter)
    if centered:
        xx -= xcenter
        yy -= ycenter
    else:
        xx += x0
        yy += y0
    geomdata["x"] = xx[boolrefs]
    geomdata["y"] = yy[boolrefs]
    geomdata["angle"] = angle[boolrefs]
    return geomdata, radius
