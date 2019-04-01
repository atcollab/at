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
import itertools
from warnings import warn
from fnmatch import fnmatch
from at.lattice import elements


class AtError(Exception):
    pass


class AtWarning(UserWarning):
    pass


def uint32_refpts(refpts, n_elements):
    """
    Return a uint32 numpy array with contents as the indices of the selected
    elements.  This is used for indexing a lattice using explicit indices.
    """
    refs = numpy.asarray(refpts).reshape(-1)
    if (refpts is None) or (refs.size is 0):
        return numpy.array([], dtype=numpy.uint32)
    elif (refs.size > n_elements+1):
        raise ValueError('too many reftps given')
    elif numpy.issubdtype(refs.dtype, numpy.bool_):
        return numpy.flatnonzero(refs).astype(numpy.uint32)

    # Handle negative indices
    refs = numpy.array([i if (i == n_elements) else i % n_elements for i in refs],
                       dtype=numpy.uint32)

    # Check ascending
    if refs.size > 1:
        prev = refs[0]
        for next in refs[1:]:
            if next < prev:
                raise ValueError('refpts should be given in ascending order')
            elif next == prev:
                raise ValueError('refpts contains duplicates or index(es) out'
                                 ' of range')
            prev = next

    return refs


def bool_refpts(refpts, n_elements):
    """
    Return a boolean numpy array of length n_elements + 1 where True elements
    are selected. This is used for indexing a lattice using True or False
    values.
    """
    if isinstance(refpts, numpy.ndarray) and refpts.dtype == bool:
        diff = 1 + n_elements - refpts.size
        if diff is 0:
            return refpts
        elif diff > 0:
            return numpy.append(refpts, numpy.zeros(diff, dtype=bool))
        else:
            return refpts[:n_elements+1]
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
    """
    if isinstance(key, elements.Element):
        elems = [elem for elem in ring if isinstance(elem, type(key))]
    elif isinstance(key, type):
        elems = [elem for elem in ring if isinstance(elem, key)]
    elif numpy.issubdtype(type(key), numpy.str_):
        elems = [elem for elem in ring if fnmatch(elem.FamName, key)]
        if not quiet:
            matched_fams = set(elem.FamName for elem in elems)
            ending = 'y' if len(matched_fams) == 1 else 'ies'
            print("String '{0}' matched {1} famil{2}: {3}\n"
                  "all corresponding elements have been "
                  "returned.".format(key, len(matched_fams), ending,
                                     ', '.join(matched_fams)))
    else:
        raise TypeError("Invalid key type {0}; please enter a string, element"
                        " type, or element instance.".format(type(key)))
    return elems


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


def get_ring_energy(ring):
    """Establish the energy of the ring from the Energy attribute of the
    elements. Energies of RingParam elements are most prioritised, if none are
    found then the energies from RFCavity elements will be used, if none are
    found then the energies from all elements will be used. An error will be
    raised if no elements have a 'Energy' attribute or if inconsistent values
    for energy are found.

    Args:
        ring: sequence of elements of which you wish to establish the energy.
    """
    rp_energies = []
    rf_energies = []
    energies = []
    for elem in ring:
        if hasattr(elem, 'Energy'):
            energies.append(elem.Energy)
            if isinstance(elem, elements.RingParam):
                rp_energies.append(elem.Energy)
            elif isinstance(elem, elements.RFCavity):
                rf_energies.append(elem.Energy)
    if not energies:
        raise AtError('Lattice energy is not defined.')
    elif rp_energies:
        energy = max(rp_energies)
    elif rf_energies:
        energy = max(rf_energies)
    else:
        energy = max(energies)
    if len(set(energies)) > 1:
        warn(AtWarning('Inconsistent energy values in ring, {0} has been '
                       'used.'.format(energy)))
    return energy


def tilt_elem(elem, rots):
    """
    set a new tilt angle to an element. The rotation matrices are stored in the R1 and R2 attributes
    R1 = [[ cos(rots) sin(rots)]    R2 = [[cos(rots) -sin(rots)]
          [-sin(rots) cos(rots)]]         [sin(rots)  cos(rots)]]
    
    :param elem: element to be tilted
    :param rots: tilt angle (in radians).
           rots > 0 corresponds to a corkskew rotation of the element looking in the direction of the beam
    :return: None
    """
    cs = numpy.cos(rots)
    sn = numpy.sin(rots)
    rm = numpy.diag([cs, cs, cs, cs, 1.0, 1.0]).T  # transpose to get Fortran alignment
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
