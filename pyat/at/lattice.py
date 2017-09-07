"""
Helper functions for working with AT lattices.

A lattice as understood by pyAT is any sequence of elements.  These functions
are useful for working with these sequences.
"""
import numpy
import collections


def uint32_refpts(refpts, n_elements):
    if isinstance(refpts, numpy.ndarray) and refpts.dtype == bool:
        urefpts = numpy.array(numpy.flatnonzero(refpts), dtype=numpy.uint32)
    else:
        if not isinstance(refpts, (collections.Sequence, numpy.ndarray)):
            refpts = [refpts]
        if numpy.any(numpy.diff(numpy.array(refpts)) < 0) or (refpts[-1] > n_elements):
            raise ValueError('refpts must be ascending and less or equal to {}'.format(n_elements))
        urefpts = numpy.asarray(refpts, dtype=numpy.uint32)
    return urefpts


def bool_refpts(refpts, n_elements):
    if isinstance(refpts, numpy.ndarray) and refpts.dtype == bool:
        return refpts
    else:
        brefpts = numpy.zeros(n_elements+1, dtype=bool)
        brefpts[refpts] = True
        return brefpts


def get_s_pos(ring, refpts):
    """
    Return a numpy array corresponding to the s position of the specified
    elements.
    """
    s_pos = numpy.concatenate(([0.0], numpy.cumsum([getattr(el, 'Length', 0.0) for el in ring])))
    return s_pos[refpts]
