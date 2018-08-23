"""
Helper functions for working with AT lattices.

A lattice as understood by pyAT is any sequence of elements.  These functions
are useful for working with these sequences.

The refpts functions allow selecting a number of points in a lattice.
Indexing runs from zero (the start of the first element) to n_elements (the end of the last element).
"""
import numpy
import itertools


def uint32_refpts(refpts, n_elements):
    """
    Return a uint32 numpy array with contents as the indices of the selected
    elements.  This is used for indexing a lattice using explicit indices.
    """
    if isinstance(refpts, numpy.ndarray) and refpts.dtype == bool:
        refs = numpy.flatnonzero(refpts)
    else:
        refs = numpy.ravel(numpy.asarray(refpts))

    if refs.size > 0:
        if (numpy.any(numpy.diff(refs) < 0) or
                (refs[-1] > n_elements) or
                (refs[0]) < 0):
            error_msg = 'refpts must be ascending and less than or equal to {}'
            raise ValueError(error_msg.format(n_elements))

    return numpy.asarray(refs, dtype=numpy.uint32)


def bool_refpts(refpts, n_elements):
    """
    Return a boolean numpy array of length n_elements + 1 where True elements
    are selected. This is used for indexing a lattice using True or False
    values.
    """
    if isinstance(refpts, numpy.ndarray) and refpts.dtype == bool:
        return refpts
    else:
        brefpts = numpy.zeros(n_elements + 1, dtype=bool)
        brefpts[refpts] = True
        return brefpts


def checkattr(*args):
    """Return a function to be used as a filter. Check the presence or the value of an attribute. This function can be
    used to extract from a ring all elements have a given attribute.

    filtfunc = checkattr(attrname)              returns the function filtfunc such that
        ok=filtfunc(element) is true if the element has a 'attrname' attribute

    filtfunc = checkattr(attrname, attrvalue)   returns the function filtfunc such that
        ok=filtfunc(element) is true if the element has a 'attrname' attribute with the value attrvalue

    Examples:

    cavs = filter(checkattr('Frequency'), ring)
        returns an iterator over all elements in ring that have a 'Frequency' attribute

    elts = filter(checkattr('K', 0.0), ring)
        returns an iterator over all elements in ring that have a 'K' attribute equal to 0.0

    """
    def testf(el):
        try:
            v = getattr(el, args[0])
            return len(args) == 1 or v == args[1]
        except AttributeError:
            return False
    return testf


def checktype(eltype):
    """Return a function to be used as a filter. Check the type of an element. This function can be
    used to extract from a ring all elements have a given type.

    filtfunc = checktype(class)                 returns the function filtfunc such that
        ok=filtfunc(element) is true if the element is an instance of class

    Example:

    qps = filter(checktype(at.Quadrupole), ring) returns an iterator over all quadrupoles in ring

    """
    return lambda el: isinstance(el, eltype)


def get_cells(ring, *args):
    """Return a numpy array of booleans, with the same length as ring, marking all elements satisfying a given condition.

    refpts = getcells(ring, filtfunc) selects all elements for which the function filtfunc(element) returns True

    refpts = getcells(ring, attrname) selects all elements having a 'attrname' attribute

    refpts = getcells(ring, attrname, attrvalue) selects all elements having a 'attrname' attribute with value attrvalue

    Example:

    refpts = getcells(ring, 'Frequency')
        returns a numpy array of booleans where all elements having a 'Frequency' attribute are True

    refpts = getcells(ring, 'K', 0.0)
        returns a numpy array of booleans where all elements having a 'K' attribute equal to 0.0 are True

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
    2) a numpy array of booleans as long as ring where selected elements are true

    """
    if isinstance(refpts, numpy.ndarray) and refpts.dtype == bool:
        for el, tst in zip(ring, refpts):
            if tst:
                yield el
    else:
        for i in refpts:
            yield ring[i]


def refpts_subset(ring, refpts):
    """Return a list of all elements in ring identified by refpts.

    refpts may be:

    1) a list of integers (0 indicating the first element)
    2) a numpy array of booleans as long as ring where selected elements are true

    """
    if isinstance(refpts, numpy.ndarray) and refpts.dtype == bool:
        return [el for el, tst in zip(ring, refpts) if tst]
    else:
        return [ring[i] for i in refpts]


def get_s_pos(ring, refpts=None):
    """
    Return a numpy array corresponding to the s position of the specified
    elements.
    """
    if refpts is None:
        refpts = len(ring)
    # Positions at the end of each element.
    s_pos = numpy.cumsum([getattr(el, 'Length', 0.0) for el in ring])
    # Prepend position at the start of the first element.
    s_pos = numpy.concatenate(([0.0], s_pos))
    refpts = uint32_refpts(refpts, len(ring))
    return numpy.squeeze(s_pos[refpts])


def tilt_elem(elem, rots):
    """
    set a new tilt angle to an element. The rotation matrices are stored in the R1 and R2 attributes
    R1 = [[ cos(rots) sin(rots)]    R2 = [[cos(rots) -sin(rots)]
          [-sin(rots) cos(rots)]]         [sin(rots)  cos(rots)]]
    
    :param elem: element to be tilted
    :param rots: tilt angle.
           rots > 0 corresponds to a corkskew rotation of the element looking in the direction of the beam
    :return: None
    """
    cs = numpy.cos(rots)
    sn = numpy.sin(rots)
    rm = numpy.diag([cs, cs, cs, cs, 1.0, 1.0]).T     # transpose to get Fortran alignment
    rm[0, 2] = sn
    rm[1, 3] = sn
    rm[2, 0] = -sn
    rm[3, 1] = -sn
    elem.R1 = numpy.asfortranarray(rm)
    elem.R2 = numpy.asfortranarray(rm.T)


def shift_elem(elem, deltax=0.0, deltaz=0.0):
    """
    set a new displacement vector to an element. Translation vectors are stored in the T1 and T2 attributes

    :param elem:  element to be displaced
    :param deltax: horizontal displacement of the element
    :param deltaz:  vertical displacement of the element
    :return: None
    """
    tr = numpy.array([deltax, 0.0, deltaz, 0.0, 0.0, 0.0])
    elem.T1 = -tr
    elem.T2 = tr


def set_tilt(ring, tilts):
    """Set tilts to a list of elements.

    ring:   any iterable over elements to be tilted
    tilts:  any iterable as long as ring containing tilt values or
            scalar tilt value to be applied to all elements

    """
    try:
        it = iter(tilts)
    except TypeError:
        it = itertools.repeat(tilts)
    for el, tilt in zip(ring, it):
        tilt_elem(el, tilt)


def set_shift(ring, dxs, dzs):
    """Set translations to a list of elements.

    ring:   any iterable over elements to be shifted
    dxs:    any iterable as long as ring containing horizontal displacement values or
            scalar horizontal displacement value to be applied to all elements
    dzs:    any iterable as long as ring containing vertical displacement values or
            scalar vertical displacement value to be applied to all elements

    """
    try:
        itx = iter(dxs)
    except TypeError:
        itx = itertools.repeat(dxs)
    try:
        itz = iter(dzs)
    except TypeError:
        itz = itertools.repeat(dzs)
    for el, dx, dy in zip(ring, itx, itz):
        shift_elem(el, dx, dy)
