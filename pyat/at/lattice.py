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
