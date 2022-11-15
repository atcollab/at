"""
Helper functions for working with AT lattices.

A lattice as understood by pyAT is any sequence of elements.  These functions
are useful for working with these sequences.

The ``refpts`` argument allow functions to select points in the lattice,
returned values are given at the entrance of each element specified in refpts;
``refpts`` can be:

#. an integer in the range [-len(ring), len(ring)-1] selecting the element
   according to python indexing rules. As a special case, len(ring) is
   allowed and refers to the end of the last element,
#. an ordered sequence of such integers without duplicates,
#. a sequence of booleans of maximum length len(ring)+1, where selected
   elements are True.
"""
import numpy
import functools
from typing import Callable, Optional, Sequence, Iterator, TypeVar
from typing import Union, Tuple, List
from itertools import compress
from fnmatch import fnmatch
from .elements import Element, Dipole

Refpts = TypeVar("Refpts", int, Sequence[int], bool, Sequence[bool])
ElementFilter = Callable[[Element], bool]
BoolRefpts = numpy.ndarray
Uint32Refpts = numpy.ndarray
Key = Union[type, Element, str]


__all__ = ['AtError', 'AtWarning', 'check_radiation', 'check_6d',
           'set_radiation', 'set_6d',
           'make_copy', 'uint32_refpts', 'bool_refpts',
           'checkattr', 'checktype', 'checkname',
           'get_cells', 'get_elements', 'get_refpts', 'get_s_pos',
           'refpts_count', 'refpts_len', 'refpts_iterator',
           'set_shift', 'set_tilt', 'tilt_elem', 'shift_elem',
           'get_value_refpts', 'set_value_refpts', 'Refpts',
           'get_geometry']


class AtError(Exception):
    pass


class AtWarning(UserWarning):
    pass


def check_radiation(rad: bool) -> Callable:
    """Deprecated decorator for optics functions (see :py:func:`check_6d`).

    :meta private:
    """
    return check_6d(rad)


def set_radiation(rad: bool) -> Callable:
    """Deprecated decorator for optics functions (see :py:func:`set_6d`).

    :meta private:
    """
    return set_6d(rad)


def check_6d(is_6d: bool) -> Callable:
    """Decorator for optics functions

    Wraps a function like ``func(ring, *args, **kwargs)`` where ring is a
    :py:class:`.Lattice` object. Raise :py:class:`.AtError` if ``ring.is_6d``
    is not the desired ``is_6d`` state. No test is performed if ring is not a
    :py:class:`.Lattice`.

    Parameters:
        is_6d: Desired 6D state

    Raises:
        AtError: if ``ring.is_6d`` is not ``is_6d``

    Example:

        .. code-block:: python

            @check_6d(False)
            def find_orbit4(ring, dct=None, guess=None, **kwargs):
                ...

        Raises :py:class:`.AtError` if ``ring.is_6d`` is :py:obj:`True`

    See Also:
        :py:func:`set_6d`
    """
    def radiation_decorator(func):
        @functools.wraps(func)
        def wrapper(ring, *args, **kwargs):
            ringrad = getattr(ring, 'is_6d', is_6d)
            if ringrad != is_6d:
                raise AtError('{0} needs "ring.is_6d" {1}'.format(
                    func.__name__, is_6d))
            return func(ring, *args, **kwargs)
        return wrapper
    return radiation_decorator


def set_6d(is_6d: bool) -> Callable:
    """Decorator for optics functions

    Wraps a function like ``func(ring, *args, **kwargs)`` where ``ring`` is a
    :py:class:`.Lattice` object. If ``ring.is_6d`` is not the desired ``is_6d``
    state, ``func`` will be called with a modified copy of ``ring`` satisfying
    the ``is_6d`` state.

    Parameters:
        is_6d: Desired 6D state

    Example:

        .. code-block:: python

            @set_6d(True)
            def compute(ring)
                ...

        is roughly equivalent to:

        .. code-block:: python

            if ring.is_6d:
                compute(ring)
            else:
                compute(ring.enable_6d(copy=True))

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
    """Decorator for optics functions

    Wraps a function like ``func(ring, refpts, *args, **kwargs)`` where
    ``ring`` is a :py:class:`.Lattice` object.

    If ``copy`` is :py:obj:`False`, the function is not modified,

    If ``copy`` is :py:obj:`True`, a shallow copy of ``ring`` is done, then the
    elements selected by ``refpts`` are deep-copied, then ``func`` is applied
    to the copy, and the new ring is returned.

    Parameters:
        copy:   If :py:obj:`True`, apply the decorated function to a copy of
          ``ring``
    """
    if copy:
        def copy_decorator(func):
            @functools.wraps(func)
            def wrapper(ring, refpts, *args, **kwargs):
                ring = ring.replace(refpts)
                func(ring, refpts, *args, **kwargs)
                return ring
            return wrapper
    else:
        def copy_decorator(func):
            return func
    return copy_decorator


def uint32_refpts(refpts: Refpts, n_elements: int) -> Uint32Refpts:
    """Returns an integer array of element indices, selecting ring elements.

    Parameters:
        refpts:     refpts may be:

          #. an integer or a sequence of integers
             (0 indicating the first element)
          #. a sequence of booleans marking the selected elements
        n_elements: Length of the lattice

    Returns:
        uint32_ref (Uint32Refpts): uint32 numpy array used for indexing
          ``Elements`` in a lattice.
    """
    refs = numpy.ravel(refpts)
    if (refpts is None) or (refs.size == 0):
        return numpy.array([], dtype=numpy.uint32)
    elif numpy.issubdtype(refs.dtype, numpy.bool_):
        return numpy.flatnonzero(refs).astype(numpy.uint32)

    # Handle negative indices
    refs = numpy.array([i if (i == n_elements) else i % n_elements
                        for i in refs], dtype=numpy.uint32)

    # Check ascending
    if refs.size > 1:
        prev = refs[0]
        for nxt in refs[1:]:
            if nxt < prev:
                raise ValueError('refpts should be given in ascending order')
            elif nxt == prev:
                raise ValueError('refpts contains duplicates or index(es) out'
                                 ' of range')
            prev = nxt

    return refs


# Private function accepting a callable for refpts
def _uint32_refs(ring: Sequence[Element], refpts: Refpts) -> Uint32Refpts:
    """Returns an integer array of element indices, selecting ring elements.

    Parameters:
        refpts:     refpts may be:

          #. an integer or a sequence of integers
             (0 indicating the first element)
          #. a sequence of booleans marking the selected elements
          #. a callable ``filtfunc`` such that ``filtfunc(elem)`` is True for
             selected elements

    Returns:
        uint32_ref (Uint32Refpts): uint32 numpy array used for indexing
          ``Elements`` in a lattice.
    """
    if callable(refpts):
        return numpy.array([i for i, el in enumerate(ring) if refpts(el)],
                           dtype=numpy.uint32)
    elif refpts is None:
        return numpy.array([], dtype=numpy.uint32)
    else:
        return uint32_refpts(refpts, len(ring))


def bool_refpts(refpts: Refpts, n_elements: int) -> BoolRefpts:
    """Returns n bool array of element indices, selecting ring elements.

    Parameters:
        refpts:     refpts may be:

          #. an integer or a sequence of integers
             (0 indicating the first element)
          #. a sequence of booleans marking the selected elements
        n_elements: Length of the lattice

    Returns:
        bool_refs (BoolRefpts):  A bool numpy array used for indexing
          ``Elements`` in a lattice.
    """
    refs = numpy.ravel(refpts)
    if (refpts is None) or (refs.size == 0):
        return numpy.zeros(n_elements + 1, dtype=bool)
    elif refs.dtype == bool:
        diff = 1 + n_elements - refs.size
        if diff <= 0:
            return refs[:n_elements + 1]
        else:
            return numpy.append(refs, numpy.zeros(diff, dtype=bool))
    else:
        brefpts = numpy.zeros(n_elements + 1, dtype=bool)
        brefpts[refs] = True
        return brefpts


# Private function accepting a callable for refpts
def _bool_refs(ring: Sequence[Element], refpts: Refpts) -> BoolRefpts:
    """Returns a bool array of element indices, selecting ring elements.

    Parameters:
        refpts:     refpts may be:

          #. an integer or a sequence of integers
             (0 indicating the first element)
          #. a sequence of booleans marking the selected elements
          #. a callable ``filtfunc`` such that ``filtfunc(elem)`` is True for
             selected elements

    Returns:
        bool_refs (BoolRefpts):  A bool numpy array used for indexing
          ``Elements`` in a lattice.
    """
    if callable(refpts):
        return numpy.array([refpts(el) for el in ring] + [False], dtype=bool)
    elif refpts is None:
        return numpy.zeros(len(ring) + 1, dtype=bool)
    else:
        return bool_refpts(refpts, len(ring))


def checkattr(attrname: str, attrvalue: Optional = None) \
        -> ElementFilter:
    # noinspection PyUnresolvedReferences
    """Checks the presence or the value of an attribute

        Returns a function to be used as an ``Element`` filter, which checks the
        presence or the value of an attribute of the provided ``Element``.
        This function can be used to extract from a ring all elements
        having a given attribute.

        Parameters:
            attrname: Attribute name
            attrvalue: Attribute value. If absent, the returned function checks
              the presence of an ``attrname`` attribute. If present, the
              returned function checks if ``attrname`` == ``attrvalue``.

        Returns:
            checkfun (ElementFilter):   Element filter function

        Examples:

            >>> cavs = filter(checkattr('Frequency'), ring)

            Returns an iterator over all elements in ring that have a
            ``Frequency`` attribute

            >>> elts = filter(checkattr('K', 0.0), ring)

            Returns an iterator over all elements in ring that have a ``K``
            attribute equal to 0.0
        """

    def testf(el):
        try:
            v = getattr(el, attrname)
            return (attrvalue is None) or (v == attrvalue)
        except AttributeError:
            return False

    return testf


def checktype(eltype: Union[type, Tuple[type, ...]]) -> ElementFilter:
    # noinspection PyUnresolvedReferences
    """Checks the type of an element

    Returns a function to be used as an ``Element`` filter, which checks the
    type of the provided ``Element``.
    This function can be used to extract from a ring all elements
    having a given type.

    Parameters:
        eltype: Desired ``Element`` type

    Returns:
        checkfun (ElementFilter):   Element filter function

    Examples:

        >>> qps = filter(checktype(at.Quadrupole), ring)

        Returns an iterator over all quadrupoles in ring
    """
    return lambda el: isinstance(el, eltype)


def checkname(pattern: str) -> ElementFilter:
    # noinspection PyUnresolvedReferences
    """Checks the name of an element

        Returns a function to be used as an ``Element`` filter, which checks the
        name of the provided ``Element``.
        This function can be used to extract from a ring all elements
        having a given name.

        Parameters:
            pattern: Desired ``Element`` name. Unix shell-style wildcards are
              supported (see ``fnmatch.fnmatch``)

        Returns:
            checkfun (ElementFilter):   Element filter function

        Examples:

            >>> qps = filter(checkname('QF*'), ring)

            Returns an iterator over all with name starting with ``QF``.
        """
    return lambda el: fnmatch(el.FamName, pattern)


# noinspection PyIncorrectDocstring
def get_cells(ring: Sequence[Element], *args) -> BoolRefpts:
    # noinspection PyUnresolvedReferences
    """
    get_cells(ring, filtfunc)
    get_cells(ring, attrname)
    get_cells(ring, attrname, attrvalue)
    Returns a bool array of element indices, selecting ring elements.

    Parameters:
        ring (Sequence[Element]):       Lattice description
        filtfunc (ElementFilter):   Filter function. Selects ``Elements``
          satisfying the filter function
        attrname (str):   Attribute name
        attrvalue (Any):  Attribute value. If absent, select the
          presence of an ``attrname`` attribute. If present, select
          ``Elements`` with ``attrname`` == ``attrvalue``.

    Returns:
        bool_refs (BoolRefpts):  numpy Array of ``bool`` with length
          len(ring)+1

    Examples:

        >>> refpts = get_cells(ring, 'Frequency')

        Returns a numpy array of booleans where all elements having a
        ``Frequency`` attribute are True

        >>> refpts = get_cells(ring, 'K', 0.0)

        Returns a numpy array of booleans where all elements having a ``K``
        attribute equal to 0.0 are True
    """
    testf = args[0] if callable(args[0]) else checkattr(*args)
    return numpy.append(numpy.fromiter(map(testf, ring), dtype=bool,
                                       count=len(ring)), False)


def refpts_iterator(ring: Sequence[Element], refpts: Refpts) \
        -> Iterator[Element]:
    """Return an iterator over selected elements in a lattice

    Parameters:
        ring:       Lattice description
        refpts:     refpts may be:

          #. an integer or a sequence of integers
             (0 indicating the first element)
          #. a sequence of booleans marking the selected elements
          #. a callable ``filtfunc`` such that ``filtfunc(elem)`` is True for
             selected elements

    Returns:
        elem_iter (Iterator[Element]):  Iterator over the elements in ``ring``
          selected by ``refpts``.
    """
    if refpts is None:
        return iter(())
    elif callable(refpts):
        return filter(refpts, ring)
    else:
        refs = numpy.ravel(refpts)
        if refs.size == 0:
            return iter(())
        elif refs.dtype == bool:
            return compress(ring, refs)
        else:
            return (ring[i] for i in refs)


# noinspection PyUnusedLocal
def refpts_count(refpts: Refpts, n_elements: int) -> int:
    """Returns the number of reference points

    Parameters:
        refpts:     refpts may be:

          #. an integer or a sequence of integers
             (0 indicating the first element)
          #. a sequence of booleans marking the selected elements
        n_elements: Lattice length

    Returns:
        nrefs (int):  The number of reference points
    """
    refs = numpy.ravel(refpts)
    if (refpts is None) or (refs.size == 0):
        return 0
    elif refs.dtype == bool:
        return numpy.count_nonzero(refs)
    else:
        return len(refs)


def refpts_len(ring: Sequence[Element], refpts: Refpts) -> int:
    """Returns the number of reference points

    Parameters:
        ring:       Lattice description
        refpts:     refpts may be:

          #. an integer or a sequence of integers
             (0 indicating the first element)
          #. a sequence of booleans marking the selected elements
          #. a callable ``filtfunc`` such that ``filtfunc(elem)`` is True for
             selected elements

    Returns:
        nrefs (int):  The number of reference points
    """
    if refpts is None:
        return 0
    elif callable(refpts):
        return len(list(filter(refpts, ring)))
    else:
        return refpts_count(refpts, len(ring))


def get_refpts(ring: Sequence[Element], key: Key, quiet=True) -> Uint32Refpts:
    """Returns an ``int`` array of element indices, selecting ring elements.

    Parameters:
        ring:   Lattice description
        key:    Element selection key. May be:

                  #. an element instance, will return all elements of the same
                     type in the lattice, e.g. key=Drift('d1', 1.0)
                  #. an element type, will return all elements of that type in
                     the lattice, e.g. ``at.elements.Sextupole``
                  #. a string to match against elements' ``FamName``, supports
                     Unix shell-style wildcards, e.g. ``'BPM_*1'``
        quiet:  if false print information about matched elements for FamName
               matches.

    Returns:
        uint32_refs (Uint32Refs):    numpy Array of ``int`` as long as the
          number of refpts

    See also:
        ``get_cells``
    """
    if isinstance(key, Element):
        checkfun = checktype(type(key))
    elif isinstance(key, type):
        checkfun = checktype(key)
    elif numpy.issubdtype(type(key), numpy.str_):
        checkfun = checkname(key)
        if not quiet:
            matched_fams = set(elem.FamName for elem in filter(checkfun, ring))
            ending = 'y' if len(matched_fams) == 1 else 'ies'
            print("String '{0}' matched {1} famil{2}: {3}\n"
                  "all corresponding elements have been "
                  "returned.".format(key, len(matched_fams), ending,
                                     ', '.join(matched_fams)))
    else:
        raise TypeError("Invalid key type {0}; please enter a string, element"
                        " type, or element instance.".format(type(key)))
    return uint32_refpts(list(map(checkfun, ring)), len(ring))


def get_elements(ring: Sequence[Element], key: Key, quiet=True) \
        -> list:
    """Returns a list of elements selected by ``key``.

    Parameters:
        ring:   Lattice description
        key:    Element selection key. May be:

                  #. an element instance, will return all elements of the same
                     type in the lattice, e.g. key=Drift('d1', 1.0)
                  #. an element type, will return all elements of that type in
                     the lattice, e.g. ``at.elements.Sextupole``
                  #. a string to match against elements' ``FamName``, supports
                     Unix shell-style wildcards, e.g. ``'BPM_*1'``
        quiet: if false print information about matched elements for FamName
               matches.

    Returns:
        elem_list (list):  list of ``Elements`` matching key
    """
    if isinstance(key, Element):
        checkfun = checktype(type(key))
    elif isinstance(key, type):
        checkfun = checktype(key)
    elif numpy.issubdtype(type(key), numpy.str_):
        checkfun = checkname(key)
        if not quiet:
            matched_fams = set(elem.FamName for elem in filter(checkfun, ring))
            ending = 'y' if len(matched_fams) == 1 else 'ies'
            print("String '{0}' matched {1} famil{2}: {3}\n"
                  "all corresponding elements have been "
                  "returned.".format(key, len(matched_fams), ending,
                                     ', '.join(matched_fams)))
    else:
        raise TypeError("Invalid key type {0}; please enter a string, element"
                        " type, or element instance.".format(type(key)))
    return list(filter(checkfun, ring))


def get_value_refpts(ring: Sequence[Element], refpts: Refpts,
                     attrname: str, index: Optional[int] = None):
    """Extracts attribute values from selected lattice ``Elements``.

    Parameters:
        ring:       Lattice description
        refpts:     refpts may be:

          #. an integer or a sequence of integers
             (0 indicating the first element)
          #. a sequence of booleans marking the selected elements
          #. a callable ``filtfunc`` such that ``filtfunc(elem)`` is True for
             selected elements
        attrname:   Attribute name
        index:      index of the value to retrieve if ``attrname`` is
          an array.

          If ``None`` the full array is retrieved

    Returns:
        attrvalues: numpy Array of attribute values.
    """
    if index is None:
        def getf(elem):
            return getattr(elem, attrname)
    else:
        def getf(elem):
            return getattr(elem, attrname)[index]

    return numpy.array([getf(elem) for elem in refpts_iterator(ring, refpts)])


def set_value_refpts(ring: Sequence[Element], refpts: Refpts,
                     attrname: str, attrvalues, index: Optional[int] = None,
                     increment: Optional[bool] = False,
                     copy: Optional[bool] = False):
    """Set the values of an attribute of an array of elements based on
    their refpts


    Parameters:
        ring:       Lattice description
        refpts:     refpts may be:

          #. an integer or a sequence of integers
             (0 indicating the first element)
          #. a sequence of booleans marking the selected elements
          #. a callable ``filtfunc`` such that ``filtfunc(elem)`` is True for
             selected elements
        attrname:   Attribute name
        attrvalues: Attribute values
        index:      index of the value to set if ``attrname`` is
          an array. if ``None``, the full array is replaced by ``attrvalue``
        increment:  Add values to the initial values.

          If False the initial value is replaced (Default)
        copy:       If ``False``, the modification is done in-place,

          If ``True``, return a shallow copy of the lattice. Only the
          modified elements are copied.

          .. Caution::

             a shallow copy means that all non-modified
             elements are shared with the original lattice.
             Any further modification will affect in both lattices.
    """
    if index is None:
        def setf(elem, value):
            setattr(elem, attrname, value)
    else:
        def setf(elem, value):
            getattr(elem, attrname)[index] = value

    if increment:
        attrvalues += get_value_refpts(ring, refpts, attrname, index=index)
    else:
        attrvalues = numpy.broadcast_to(attrvalues, (refpts_len(ring, refpts),))

    # noinspection PyShadowingNames
    @make_copy(copy)
    def apply(ring, refpts, values):
        for elm, val in zip(refpts_iterator(ring, refpts), values):
            setf(elm, val)

    return apply(ring, refpts, attrvalues)


def get_s_pos(ring: Sequence[Element], refpts: Optional[Refpts] = None) \
        -> Sequence[float]:
    """Returns the locations of selected elements

    Parameters:
        ring:       Lattice description
        refpts:     refpts may be:

          #. an integer or a sequence of integers
             (0 indicating the first element)
          #. a sequence of booleans marking the selected elements
          #. a callable ``filtfunc`` such that ``filtfunc(elem)`` is True for
             selected elements

    Returns:
        s_pos:  Array of locations of the elements selected by ``refpts``
    """
    if refpts is None:
        refpts = range(len(ring) + 1)
    # Positions at the end of each element.
    s_pos = numpy.cumsum([getattr(el, 'Length', 0.0) for el in ring])
    # Prepend position at the start of the first element.
    s_pos = numpy.concatenate(([0.0], s_pos))
    refpts = _uint32_refs(ring, refpts)
    return s_pos[refpts]
    
    
def rotate_elem(elem: Element, tilt: float = 0.0, pitch: float = 0.0,
                yaw: float = 0.0, relative: bool = False) -> None:
    r"""Set the tilt, pitch and yaw angle of an ``Element``.
    The tilt is a rotation around the ``s`` axis, the pitch is a
    rotation around the ``x``-axis and the yaw is a rotation around
    the ``y``-axis.
    A positive angle represent a clockwise rotation when looking in
    the direction of the roation axis.

    Parameters:
        elem:           Element to be tilted
        tilt:           Tilt angle [rd]
        pitch:          Pitch angle [rd]
        yaw:            Yaw angle [rd]
        relative:       If True, the rotation is added to the previous one
    """ 
    def _get_rm_tv(le, tilt, pitch, yaw):
        ct = numpy.cos(tilt)
        st = numpy.sin(tilt)
        ap, ay = -0.5*le*numpy.sin(pitch), -0.5*le*numpy.sin(yaw)      
        rr1 = numpy.asfortranarray(numpy.diag([ct, ct, ct, ct, 1.0, 1.0]))
        rr1[0, 2] = st
        rr1[1, 3] = st 
        rr1[2, 0] = -st
        rr1[3, 1] = -st
        rr2 = tm1.T      
        t1 = numpy.array([ay, -yaw, ap, -pitch, 0, 0])
        t2 = numpy.array([ay, yaw, ap, pitch, 0, 0])     
        rt1 = numpy.asfortranarray(numpy.diag(numpy.ones(6)))
        rt1[1, 4] =  ct*t1[1]
        rt1[3, 4] =  ct*t1[3]       
        rt2 = numpy.asfortranarray(numpy.diag(numpy.ones(6)))
        rt2[1, 4] =  ct*t2[1]
        rt2[3, 4] =  ct*t2[3]
        return rr1 @ rt1, rt2 @ rr2, t1, t2
        
    rr10 = numpy.asfortranarray(numpy.diag(numpy.ones(6)))
    rr10[:4, :4] = elem.R1[:4, :4]
    rt10 = rr10.T @ elem.R1     
    tilt0 = numpy.arctan2(rr10[0, 2], rr10[0, 0])
    yaw0 = rt10[1, 4]/rr10[0, 0]
    pitch0 = rt10[3, 4]/rr10[0, 0]
    _, _, t10, t20 = _get_rm_tv(elem.Length, tilt0, pitch0, yaw0)    
  
  
    if relative and hasattr(elem, 'R1') and hasattr(elem, 'R2'):
 
        
    r1, r2, t1, t2 = _get_rm_tv(elem.Length, tilt, pitch, yaw)      
            
    if relative and hasattr(elem, 'T1') and hasattr(elem, 'T2'):
        t1 -= elem.T1
        t2 += elem.T2       
     
    elem.R1 = r1
    elem.R2 = r2
    elem.T1 = t1
    elem.T2 = t2


def tilt_elem(elem: Element, rots: float, relative: bool = False) -> None:
    r"""Set the tilt angle :math:`\theta` of an ``Element``

    The rotation matrices are stored in the ``R1`` and ``R2`` attributes

    :math:`R_1=\begin{pmatrix} cos\theta & sin\theta \\
    -sin\theta & cos\theta \end{pmatrix}`,
    :math:`R_2=\begin{pmatrix} cos\theta & -sin\theta \\
    sin\theta & cos\theta \end{pmatrix}`

    Parameters:
        elem:           Element to be tilted
        rots:           Tilt angle :math:`\theta` [rd].
          ``rots`` > 0 corresponds to a corkscrew rotation of the element
          looking in the direction of the beam
        relative:       If True, the rotation is added to the previous one
    """
    cs = numpy.cos(rots)
    sn = numpy.sin(rots)
    rm = numpy.asfortranarray(numpy.diag([cs, cs, cs, cs, 1.0, 1.0]))
    rm[0, 2] = sn
    rm[1, 3] = sn
    rm[2, 0] = -sn
    rm[3, 1] = -sn
    if relative and hasattr(elem, 'R1') and hasattr(elem, 'R2'):
        elem.R1 = rm @ elem.R1
        elem.R2 = rm.T @ elem.R2
    else:
        elem.R1 = rm
        elem.R2 = rm.T


def shift_elem(elem: Element, deltax: float = 0.0, deltaz: float = 0.0,
               relative: Optional[bool] = False) -> None:
    """Sets the transverse displacement of an ``Element``

    The translation vectors are stored in the ``T1`` and ``T2`` attributes

    Parameters:
        elem:           Element to be shifted
        deltax:         Horizontal displacement [m]
        deltaz:         Vertical displacement [m]
        relative:       If True, the translation is added to the previous one
    """
    tr = numpy.array([deltax, 0.0, deltaz, 0.0, 0.0, 0.0])
    if relative and hasattr(elem, 'T1') and hasattr(elem, 'T2'):
        elem.T1 -= tr
        elem.T2 += tr
    else:
        elem.T1 = -tr
        elem.T2 = tr


def set_tilt(ring: Sequence[Element], tilts, relative=False) -> None:
    """Sets the tilts of a list of elements.

    Parameters:
        ring:           Lattice description
        tilts:          Sequence of tilt values as long as ring or
          scalar tilt value applied to all elements
        relative:       If True, the rotation is added to the previous one
    """
    tilts = numpy.broadcast_to(tilts, (len(ring),))
    for el, tilt in zip(ring, tilts):
        tilt_elem(el, tilt, relative=relative)


def set_shift(ring: Sequence[Element], dxs, dzs, relative=False) -> None:
    """Sets the translations of a list of elements.

    Parameters:
        ring:           Lattice description
        dxs:            Sequence of horizontal displacements values as long as
          ring or scalar value applied to all elements [m]
        dzs:            Sequence of vertical displacements values as long as
          ring or scalar value applied to all elements [m]
        relative:       If True, the displacement is added to the previous one
    """
    dxs = numpy.broadcast_to(dxs, (len(ring),))
    dzs = numpy.broadcast_to(dzs, (len(ring),))
    for el, dx, dy in zip(ring, dxs, dzs):
        shift_elem(el, dx, dy, relative=relative)


def get_geometry(ring: List[Element],
                 start_coordinates: Tuple[float, float, float] = (0, 0, 0),
                 centered: bool = False):
    """Compute the 2D ring geometry in cartesian coordinates

    Parameters:
        ring: Lattice description
        start_coordinates: x,y,angle at starting point
        centered: it True the coordinates origin is the center of the ring

    Returns:
        geomdata: recarray containing, x, y, angle
        radius: machine radius

    Example:

       >>> geomdata, radius = get_geometry(ring)
    """

    geom_dtype = [('x', numpy.float64, (1, )),
                  ('y', numpy.float64, (1, )),
                  ('angle', numpy.float64, (1, ))]
    geomdata = numpy.recarray((len(ring)+1, ), dtype=geom_dtype)
    xx = numpy.zeros((len(ring)+1, 1))
    yy = numpy.zeros((len(ring)+1, 1))
    angle = numpy.zeros((len(ring)+1, 1))
    x, y, t = start_coordinates
    x0, y0, t0 = start_coordinates

    for ind, el in enumerate(ring+[ring[0]]):
        ll = el.Length
        if isinstance(el, Dipole) and el.BendingAngle != 0:
            ang = 0.5 * el.BendingAngle
            ll *= numpy.sin(ang)/ang
        else:
            ang = 0.0
        t -= ang
        x += ll * numpy.cos(t)
        y += ll * numpy.sin(t)
        t -= ang
        xx[ind] = x
        yy[ind] = y
        angle[ind] = t

    radius = get_s_pos(ring, len(ring))[0] / abs(t-t0) \
        if t != 0.0 else 0.0
    if centered:
        xx += -abs(radius)*numpy.sin(t0) - x0
        yy += abs(radius)*numpy.cos(t0) - y0
    geomdata['x'] = xx
    geomdata['y'] = yy
    geomdata['angle'] = angle
    return geomdata, radius
