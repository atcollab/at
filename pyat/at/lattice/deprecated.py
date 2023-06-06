from .lattice_object import Lattice
from .elements import Element
from .utils import Refpts, BoolRefpts, Uint32Refpts
from .utils import checkattr, get_uint32_index, get_bool_index
from typing import Sequence

__all__ = ['get_cells', 'get_refpts']


# noinspection PyIncorrectDocstring
def get_cells(ring: Sequence[Element], refpts: Refpts, *args,
              regex=False) -> BoolRefpts:
    # noinspection PyShadowingNames
    r"""
    get_cells(ring, filtfunc) -> BoolRefpts
    get_cells(ring, element_type) -> BoolRefpts
    get_cells(ring, attrname) -> BoolRefpts
    get_cells(ring, attrname, attrvalue) -> BoolRefpts
    Returns a bool array of element indices, selecting ring elements.

    Deprecated: :pycode:`get_cells(ring, refpts)` is
    :pycode:`ring.bool_refpts(refpts)` except for :py:obj:`str` arguments:
    :pycode:`get_cells(ring, attrname [, attrvalue])` is
    :pycode:`ring.bool_refpts(checkattr(strkey [, attrvalue]))`

    Parameters:
        ring (Sequence[Element]):       Lattice description
        filtfunc (ElementFilter):       Filter function. Selects
          :py:class:`.Element`\ s satisfying the filter function
        element_type (Type[Element]):   Element type
        attrname (str):                 Attribute name
        attrvalue (Any):                Attribute value. If absent, select the
          presence of an *attrname* attribute. If present, select
          :py:class:`.Element`\ s with :pycode:`attrname == attrvalue`.
        regex: Use regular expression for refpts string matching;
            Default: False (Unix shell-style wildcards)

    Returns:
        bool_refs (BoolRefpts):  numpy Array of :py:obj:`bool` with length
          len(ring)+1

    Examples:

        >>> refpts = get_cells(ring, 'Frequency')

        Returns a numpy array of booleans where all elements having a
        :pycode:`Frequency` attribute are :py:obj:`True`

        >>> refpts = get_cells(ring, 'K', 0.0)

        Returns a numpy array of booleans where all elements having a
        :pycode:`K` attribute equal to 0.0 are :py:obj:`True`

    See also:
        :py:meth:`.Lattice.bool_refpts`, :py:meth:`.Lattice.uint32_refpts`
    """
    if isinstance(refpts, str):
        refpts = checkattr(refpts, *args)
    return get_bool_index(ring, refpts, regex=regex)


# noinspection PyUnusedLocal,PyIncorrectDocstring
def get_refpts(ring: Sequence[Element], refpts: Refpts,
               regex=False) -> Uint32Refpts:
    r"""Return a :py:obj:`~numpy.uint32` array of element indices selecting
    ring elements.

    Deprecated: :pycode:`get_elements(ring, refpts)` is
    :pycode:`ring.uint32_refpts(refpts)`

    Parameters:
        ring:           Lattice description
        refpts:         Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        regex: Use regular expression for refpts string matching;
            Default: False (Unix shell-style wildcards)

    Returns:
        uint32_refs (Uint32Refs):    :py:obj:`~numpy.uint32` numpy array as
          long as the number of refpts

    See also:
        :py:meth:`.Lattice.uint32_refpts`, :py:meth:`.Lattice.bool_refpts`
    """
    return get_uint32_index(ring, refpts, regex=regex)


Lattice.uint32_refpts = get_uint32_index
Lattice.bool_refpts = get_bool_index
Lattice.get_cells = get_cells
Lattice.get_refpts = get_refpts
