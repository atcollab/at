"""
Helper functions for working with AT lattices.

A lattice as understood by pyAT is any sequence of elements.  These functions
are useful for working with these sequences.

The refpts allow functions to select points in the lattice, returned values are
given at the entrance of each element specified in refpts;
refpts can be:
    - an integer in the range [-len(ring), len(ring)-1] selecting the element
        according to python indexing rules. As a special case, len(ring) is
        allowed and refers to the end of the last element,
    - an ordered list of such integers without duplicates,
    - a numpy array of booleans of maximum length len(ring)+1, where selected
        elements are True.
"""
import numpy
import functools
from fnmatch import fnmatch
from at.lattice import elements

__all__ = ['AtError', 'AtWarning', 'check_radiation', 'uint32_refpts',
           'bool_refpts', 'checkattr', 'checktype', 'get_cells',
           'get_elements', 'get_refpts', 'refpts_iterator', 'get_s_pos',
           'set_shift', 'set_tilt', 'tilt_elem', 'shift_elem',
           'get_value_refpts', 'set_value_refpts']


class AtError(Exception):
    pass


class AtWarning(UserWarning):
    pass


def check_radiation(rad):
    """Function to be used as a decorator for optics functions

    If ring is a Lattice object, raises an exception
        if ring.radiation is not rad

    If ring is any other sequence, no test is performed
    """
    def radiation_decorator(func):
        @functools.wraps(func)
        def wrapper(ring, *args, **kwargs):
            ringrad = getattr(ring, 'radiation', rad)
            if ringrad != rad:
                raise AtError('{0} needs radiation {1}'.format(
                    func.__name__, 'ON' if rad else 'OFF'))
            return func(ring, *args, **kwargs)
        return wrapper
    return radiation_decorator


def uint32_refpts(refpts, n_elements):
    """
    Return a uint32 numpy array with contents as the indices of the selected
    elements.  This is used for indexing a lattice using explicit indices.
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


def bool_refpts(refpts, n_elements):
    """
    Return a boolean numpy array of length n_elements + 1 where True elements
    are selected. This is used for indexing a lattice using True or False
    values.
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


def checkattr(*args):
    """Return a function to be used as a filter. Check the presence or the
    value of an attribute. This function can be used to extract from a ring
    all elements have a given attribute.

    filtfunc = checkattr(attrname)
        returns the function filtfunc such that ok=filtfunc(element) is True if
        the element has a 'attrname' attribute

    filtfunc = checkattr(attrname, attrvalue)
        returns the function filtfunc such that ok=filtfunc(element) is True if
        the element has a 'attrname' attribute with the value attrvalue

    Examples:

    cavs = filter(checkattr('Frequency'), ring)
        returns an iterator over all elements in ring that have a
        'Frequency' attribute

    elts = filter(checkattr('K', 0.0), ring)
        returns an iterator over all elements in ring that have a 'K'
        attribute equal to 0.0
    """

    def testf(el):
        try:
            v = getattr(el, args[0])
            return len(args) == 1 or v == args[1]
        except AttributeError:
            return False

    return testf


def checktype(eltype):
    """Return a function to be used as a filter. Check the type of an element.
    This function can be used to extract from a ring all elements have a
    given type.

    filtfunc = checktype(class)
        returns the function filtfunc such that ok=filtfunc(element) is True
        if the element is an instance of class

    Example:

    qps = filter(checktype(at.Quadrupole), ring)
        returns an iterator over all quadrupoles in ring
    """
    return lambda el: isinstance(el, eltype)


def checkname(pattern):
    """Return a function to be used as a filter. Check the name of an element.
    This function can be used to extract from a ring all elements have a
    given FamName attribute.

    filtfunc = checkname(pattern)
        returns the function filtfunc such that ok=filtfunc(element) is True
        if the element's FamName matches pattern.

    Example:

    qps = filter(checkname('QF.*'), ring)
        returns an iterator over all elements
    """
    return lambda el: fnmatch(el.FamName, pattern)


def get_cells(ring, *args):
    """Return a numpy array of booleans, with the same length as ring,
    marking all elements satisfying a given condition.

    refpts = getcells(ring, filtfunc)
        selects all elements for which the function filtfunc(element)
        returns True

    refpts = getcells(ring, attrname)
        selects all elements having a 'attrname' attribute

    refpts = getcells(ring, attrname, attrvalue)
        selects all elements having a 'attrname' attribute with value attrvalue

    Example:

    refpts = getcells(ring, 'Frequency')
        returns a numpy array of booleans where all elements having a
        'Frequency' attribute are True

    refpts = getcells(ring, 'K', 0.0)
        returns a numpy array of booleans where all elements having a 'K'
        attribute equal to 0.0 are True
    """
    if callable(args[0]):
        testf = args[0]
    else:
        testf = checkattr(*args)
    return numpy.array(tuple(map(testf, ring)), dtype=bool)


def refpts_iterator(ring, refpts):
    """Return an iterator over all elements in ring identified by refpts.

    refpts may be:

    1) a list of integers (0 indicating the first element)
    2) a numpy array of booleans as long as ring where selected elements
       are true
    """
    refs = numpy.ravel(refpts)
    if (refpts is None) or (refs.size == 0):
        return
    elif refs.dtype == bool:
        for el, tst in zip(ring, refs):
            if tst:
                yield el
    else:
        for i in refs:
            yield ring[i]


def refpts_count(refpts, n_elements):
    refs = numpy.ravel(refpts)
    if (refpts is None) or (refs.size == 0):
        return 0
    elif refs.dtype == bool:
        return numpy.count_nonzero(refs)
    else:
        return len(refs)


def get_refpts(ring, key, quiet=True):
    """Get the elements refpts of a family or class (type) from the lattice.

    Args:
        ring: lattice from which to retrieve the elements.
        key: can be:
             1) an element instance, will return all elements of the same type
                in the lattice, e.g. key=Drift('d1', 1.0)
             2) an element type, will return all elements of that type in the
                lattice, e.g. key=at.elements.Sextupole
             3) a string to match against elements' FamName, supports Unix
                shell-style wildcards, e.g. key='BPM_*1'
        quiet: if false print information about matched elements for FamName
               matches, defaults to True.

    Returns:
        elems: a list of elems refpts matching key
    """
    if isinstance(key, elements.Element):
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


def get_elements(ring, key, quiet=True):
    """Get the elements of a family or class (type) from the lattice.

    Args:
        ring: lattice from which to retrieve the elements.
        key: can be:
             1) an element instance, will return all elements of the same type
                in the lattice, e.g. key=Drift('d1', 1.0)
             2) an element type, will return all elements of that type in the
                lattice, e.g. key=at.elements.Sextupole
             3) a string to match against elements' FamName, supports Unix
                shell-style wildcards, e.g. key='BPM_*1'
        quiet: if false print information about matched elements for FamName
               matches, defaults to True.

    Returns:
        a list of elems matching key
    """
    if isinstance(key, elements.Element):
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


def get_value_refpts(ring, refpts, var, index=None, order=None):
    """Get the values of an attribute of an array of elements based on
    their refpts

    Args:
        ring: lattice from which to retrieve the elements.
        refpts: integer or array of integer or booleans
        var: attribute name
        index (optional): index of the value to change in case var is
        an array, if None the full array is returned (Default)
    """
    if order is not None:   # for backward compatibility with the
        index = order       # obsolete keyword 'order'

    if index is None:
        def getf(elem):
            return getattr(elem, var)
    else:
        def getf(elem):
            return getattr(elem, var)[index]

    return numpy.array([getf(elem) for elem in refpts_iterator(ring, refpts)])


def set_value_refpts(ring, refpts, var, values, index=None, order=None,
                     increment=False):
    """Set the values of an attribute of an array of elements based on
    their refpts

    Args:
        ring: lattice from which to retrieve the elements.
        refpts: integer or array of integer or booleans
        var: attribute name
        values: desired value for the attribute
        index (optional): index of the value to change in case var is an
        array, if None the full array is replaced by value (Default)
        increment (optional): adds value to the initial value, if False
        the initial value is replaced (Default)
    """
    if order is not None:   # for backward compatibility with the
        index = order       # obsolete keyword 'order'

    if index is None:
        def setf(elem, value):
            setattr(elem, var, value)
    else:
        def setf(elem, value):
            getattr(elem, var)[index] = value

    if increment:
        values = values + get_value_refpts(ring, refpts, var, index=index)
    else:
        values = values + numpy.zeros(refpts_count(refpts, len(ring)))

    for elm, val in zip(refpts_iterator(ring, refpts), values):
        setf(elm, val)


def get_s_pos(ring, refpts=None):
    """
    Return a numpy array corresponding to the s position of the specified
    elements.

    Args:
        ring: lattice from which to retrieve s position
        refpts: elements at which to return s position. If None, return
                s position at all elements in the ring.
    """
    if refpts is None:
        refpts = range(len(ring) + 1)
    # Positions at the end of each element.
    s_pos = numpy.cumsum([getattr(el, 'Length', 0.0) for el in ring])
    # Prepend position at the start of the first element.
    s_pos = numpy.concatenate(([0.0], s_pos))
    refpts = uint32_refpts(refpts, len(ring))
    return s_pos[refpts]


def tilt_elem(elem, rots, relative=False):
    """
    set a new tilt angle to an element.
    The rotation matrices are stored in the R1 and R2 attributes

    R1 = [[ cos(rots) sin(rots)]    R2 = [[cos(rots) -sin(rots)]
          [-sin(rots) cos(rots)]]         [sin(rots)  cos(rots)]]

    elem            element to be tilted
    rots            tilt angle (in radians).
                    rots > 0 corresponds to a corkskew rotation of the element
                    looking in the direction of the beam
    relative=False  Rotation relative to the previous element rotation
    """
    cs = numpy.cos(rots)
    sn = numpy.sin(rots)
    rm = numpy.asfortranarray(numpy.diag([cs, cs, cs, cs, 1.0, 1.0]))
    rm[0, 2] = sn
    rm[1, 3] = sn
    rm[2, 0] = -sn
    rm[3, 1] = -sn
    if relative and hasattr(elem, 'R1') and hasattr(elem, 'R2'):
        elem.R1 = elem.R1.dot(rm)
        elem.R2 = rm.T.dot(elem.R2)
    else:
        elem.R1 = rm
        elem.R2 = rm.T


def shift_elem(elem, deltax=0.0, deltaz=0.0, relative=False):
    """
    set a new displacement vector to an element.
    The ranslation vectors are stored in the T1 and T2 attributes

    elem            element to be displaced
    deltax          horizontal displacement of the element
    deltaz          vertical displacement of the element
    relative=False  Displacement relative to the previous alignment
    """
    tr = numpy.array([deltax, 0.0, deltaz, 0.0, 0.0, 0.0])
    if relative and hasattr(elem, 'T1') and hasattr(elem, 'T2'):
        elem.T1 -= tr
        elem.T2 += tr
    else:
        elem.T1 = -tr
        elem.T2 = tr


def set_tilt(ring, tilts, relative=False):
    """Set tilts to a list of elements.

    ring            sequence of elements to be tilted
    tilts           sequence of tilt values as long as ring or
                    scalar tilt value applied to all elements
    relative=False  Rotation relative to the previous tilt angle
    """
    tilts = numpy.broadcast_to(tilts, (len(ring),))
    for el, tilt in zip(ring, tilts):
        tilt_elem(el, tilt, relative=relative)


def set_shift(ring, dxs, dzs, relative=False):
    """Set translations to a list of elements.

    ring            sequence of elements to be shifted
    dxs             sequence of horizontal displacement as long as ring or
                    scalar horizontal displacement applied to all elements
    dzs             sequence of vertical displacement as long as ring or
                    scalar vertical displacement applied to all elements
    relative=False  Displacement relative to the previous alignment
    """
    dxs = numpy.broadcast_to(dxs, (len(ring),))
    dzs = numpy.broadcast_to(dzs, (len(ring),))
    for el, dx, dy in zip(ring, dxs, dzs):
        shift_elem(el, dx, dy, relative=relative)
