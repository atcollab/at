"""
Helper functions for working with AT lattices.

A lattice as understood by pyAT is any sequence of elements.  These functions
are useful for working with these sequences.

The refpts functions allow selecting a number of points in a lattice.
Indexing runs from zero (the start of the first element) to n_elements + 1
(the end of the last element).
"""
import numpy
import collections


def uint32_refpts(refpts, n_elements):
    """
    Return a uint32 numpy array with contents as the indices of the selected
    elements.  This is used for indexing a lattice using explicit indices.
    """
    if isinstance(refpts, numpy.ndarray) and refpts.dtype == bool:
        urefpts = numpy.array(numpy.flatnonzero(refpts), dtype=numpy.uint32)
    else:
        if not isinstance(refpts, (collections.Sequence, numpy.ndarray)):
            refpts = [refpts]
        if (numpy.any(numpy.diff(numpy.array(refpts)) < 0) or
                (refpts[-1] > n_elements) or
                numpy.any(numpy.array(refpts) < 0)):
            error_msg = 'refpts must be ascending and less than or equal to {}'
            raise ValueError(error_msg.format(n_elements))
        urefpts = numpy.asarray(refpts, dtype=numpy.uint32)
    return urefpts


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


def get_s_pos(ring, refpts=None):
    """
    Return a numpy array corresponding to the s position of the specified
    elements.
    """
    # Positions at the end of each element.
    s_pos = numpy.cumsum([getattr(el, 'Length', 0.0) for el in ring])
    # Prepend position at the start of the first element.
    s_pos = numpy.concatenate(([0.0], s_pos))
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
    :return:
    """
    tr = numpy.array([deltax, 0.0, deltaz, 0.0, 0.0, 0.0])
    elem.T1 = -tr
    elem.T2 = tr