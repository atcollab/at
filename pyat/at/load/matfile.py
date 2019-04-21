"""
Load lattices from Matlab files.
"""
from warnings import warn
import scipy.io
import numpy
from at.lattice import elements, Lattice, AtError, AtWarning
from at.load import element_from_dict

TWO_PI_ERROR = 1.E-4


def _load_element(index, element_array, check=True, quiet=False):
    """Load what scipy produces into a pyat element object.
    """
    kwargs = {}
    for field_name in element_array.dtype.fields:
        # Remove any surplus dimensions in arrays.
        data = numpy.squeeze(element_array[field_name])
        # Convert strings in arrays back to strings.
        if data.dtype.type is numpy.unicode_:
            data = str(data)
        kwargs[field_name] = data

    return element_from_dict(kwargs, index=index, check=check, quiet=quiet)


def _scanner(elems, **kwargs):
    """Extract the lattice attributes from an element list

    energy is taken from:
        1) Radiative elements (Cavity, RingParam, Radiating magnets)
        2) Any other element

    periodicity is taken from:
        1) RingParam element
        2) Sum ofthe bending angles of magnets
        3) 1

    name is taken from
        1) RingParam element
        2) ""

    All RingParam attributes are extracted
    The function removes the Energy attribute of non-radiating elements
    """
    _translate = {'Energy': 'energy', 'Periodicity': 'periodicity',
                  'FamName': 'name'}
    _ignore = {'energy', 'PassMethod', 'Length'}

    attributes = {}
    params = []
    rad_energies = []
    el_energies = []
    thetas = []

    radiate = False
    for elem in elems:
        if isinstance(elem, elements._RingParam):
            params.append(elem)
        if (isinstance(elem, (elements.RFCavity, elements._RingParam)) or
                elem.PassMethod.endswith('RadPass')):
            rad_energies.append(elem.Energy)
        elif hasattr(elem, 'Energy'):
            el_energies.append(elem.Energy)
            del elem.Energy
        if isinstance(elem, elements.Dipole):
            # noinspection PyUnresolvedReferences
            thetas.append(elem.BendingAngle)
        if (elem.PassMethod.endswith('RadPass') or
                elem.PassMethod.endswith('CavityPass')):
            radiate = True
    attributes['_radiation'] = radiate

    if 'energy' not in kwargs:
        # Guess energy from the Energy attribute of the elements
        if rad_energies:
            # Energy of radiating elements must be consistent
            energy = max(rad_energies)
            if min(rad_energies) < energy:
                raise AtError('Inconsistent energy values,')
            attributes['energy'] = energy
        elif el_energies:
            # Energy of other elements is informative only
            energy = max(el_energies)
            if min(el_energies) < energy:
                warn(AtWarning('Inconsistent energy values, '
                               '"energy" set to {0}'.format(energy)))
            attributes['energy'] = energy
        else:
            raise AtError('Lattice energy is not defined')

    if params:
        # At least one RingParam element, use the 1st one
        if len(params) > 1:
            warn(AtWarning(
                'More than 1 RingParam element, the 1st one is used'))
        attributes.update((_translate.get(key, key.lower()), attr)
                          for key, attr in vars(params[0]).items()
                          if key not in _ignore)
    else:
        # No RingParam element, try to guess
        attributes['name'] = ''
        if 'periodicity' not in kwargs:
            # Guess periodicity from the bending angles of the lattice
            try:
                nbp = 2.0 * numpy.pi / sum(thetas)
            except ZeroDivisionError:
                warn(AtWarning('No bending in the cell, '
                               'set "Periodicity" to 1'))
                attributes['periodicity'] = 1
            else:
                periodicity = int(round(nbp))
                if abs(periodicity - nbp) > TWO_PI_ERROR:
                    warn(AtWarning(
                        'Non-integer number of cells: '
                        '{0} -> {1}'.format(nbp, periodicity)))
                attributes['periodicity'] = periodicity

    return attributes


def load_mat(filename, key=None, check=True, quiet=False, keep_all=False,
             **kwargs):
    """Load a matlab at structure into a Python at list

    PARAMETERS
        filename        name of a '.mat' file
        key             name of the Matlab variable containing the lattice.
                        Default: Matlab variable name if there is only one,
                        otherwise 'RING'

    KEYWORDS
        check=True      if False, skip the coherence tests
        quiet=False     If True, suppress the warning for non-standard classes
        keep_all=False  if True, keep RingParam elements as Markers
        name            Name of the lattice
                        (default: taken from the lattice, or '')
        energy          Energy of the lattice
                        (default: taken from the elements)
        periodicity     Number of periods
                        (default: taken from the elements, or 1)

    OUTPUT
        list    pyat ring
    """

    def substitute(elem):
        if isinstance(elem, elements._RingParam):
            params = vars(elem).copy()
            name = params.pop('FamName')
            return elements.Marker(name, **params)
        else:
            return elem

    m = scipy.io.loadmat(filename)
    if key is None:
        matvars = [varname for varname in m if not varname.startswith('__')]
        key = matvars[0] if (len(matvars) == 1) else 'RING'
    element_arrays = m[key].flat
    elem_list = [_load_element(i, elem[0][0], check=check, quiet=quiet) for
                 (i, elem) in enumerate(element_arrays)]
    attrs = _scanner(elem_list, **kwargs)
    attrs.update(kwargs)
    if keep_all:
        elem_list = (substitute(elem) for elem in elem_list)
    else:
        elem_list = (elem for elem in elem_list if
                     not isinstance(elem, elements._RingParam))
    return Lattice(elem_list, **attrs)
