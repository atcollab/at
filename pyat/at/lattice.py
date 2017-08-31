"""
Helper functions for working with AT lattices.

A lattice as understood by pyAT is any sequence of elements.  These functions
are useful for working with these sequences.
"""
import numpy


def get_refpts(refpts, n_elements, append_last=False):
    """
    If using refpts in a call to at.atpass(), it must be a numpy array of type
    uint32.  This function creates such an array.
    """
    if refpts is None:
        refpts = [n_elements]
    int_array = numpy.array(refpts, dtype=numpy.uint32)
    if numpy.any(numpy.diff(int_array) < 0) or numpy.any(int_array > n_elements):
        raise ValueError('refpts must be an ascending array: {} {}'.
                         format(int_array, n_elements))
    if append_last and refpts[-1] != n_elements:
        int_array = numpy.append(int_array, [n_elements])
    return numpy.array(int_array, dtype=numpy.uint32)


def get_s_pos(ring, refpts=None):
    """
    Return a numpy array corresponding to the s position of the specified
    elements.
    """
    refpts = get_refpts(refpts, len(ring), append_last=True)
    total = 0
    s_pos = numpy.zeros(len(refpts))
    j = 1
    for i, element in enumerate(ring):
        total += element.Length
        if i in refpts:
            s_pos[j] = total
            j += 1
    s_pos[-1] = total
    return s_pos
