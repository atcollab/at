"""
Load lattice from Matlab file.

This is working from a specific file and may not be general.
"""
import scipy.io
from . import elements
import numpy


CLASS_MAPPING = {'Quad': 'Quadrupole', 'Sext': 'Sextupole', 'AP': 'Aperture',
                 'RF': 'RFCavity', 'BPM': 'Monitor'}

CLASSES = set(['Marker', 'Monitor', 'Aperture', 'Drift', 'ThinMultipole',
               'Multipole', 'Dipole', 'Bend', 'Quadrupole', 'Sextupole',
               'Octupole', 'RFCavity', 'RingParam', 'M66', 'Corrector'])

FAMILY_MAPPING = {'ap': 'Aperture', 'rf': 'RFCavity', 'bpm': 'Monitor'}

PASSMETHOD_MAPPING = {'CorrectorPass': 'Corrector', 'Matrix66Pass': 'M66',
                      'CavityPass': 'RFCavity', 'QuadLinearPass': 'Quadrupole',
                      'BendLinearPass': 'Bend', 'AperturePass': 'Aperture',
                      'ThinMPolePass': 'ThinMultipole', 'DriftPass': 'Drift'}


def hasattrs(kwargs, *attributes):
    """Check if the element would have the specified attribute(s), i.e. if they
        are in kwargs; allows checking for multiple attributes in one go.

    Args:
        kwargs (dict): The dictionary of keyword arguments passed to the Element
                        constructor.
        attributes (iterable): A list of strings, the attribute names to be
                                checked.

    Returns:
        bool: A single boolean, True if the element has any of the specified
               attributes.
    """
    results = []
    for attribute in attributes:
        try:
            kwargs[attribute]
            results.append(True)
        except KeyError:
            results.append(False)
    return any(results)


def find_class_name(kwargs):
    """Attempts to correctly identify the Class of the element from its kwargs.

    Args:
        kwargs (dict): The dictionary of keyword arguments passed to the Element
                        constructor.

    Returns:
        str: The guessed Class name, as a string.

    Raises:
        AttributeError: if the Class name found in kwargs is not a valid Class
                         in pyAT.
    """
    try:
        class_name = kwargs.pop('Class')
        class_name = CLASS_MAPPING.get(class_name, class_name)
        if class_name in CLASSES:
            return class_name
        else:
            raise AttributeError("Class {0} does not exist.".format(class_name))
    except KeyError:
        fam_name = kwargs.get('FamName')
        if fam_name in CLASSES:
            return fam_name
        elif fam_name.lower() in FAMILY_MAPPING.keys():
            return FAMILY_MAPPING[fam_name.lower()]
        else:
            pass_method = kwargs.get('PassMethod')
            class_from_pass = PASSMETHOD_MAPPING.get(pass_method)
            if class_from_pass != None:
                return class_from_pass
            else:
                if hasattrs(kwargs, 'FullGap', 'FringeInt1', 'FringeInt2', 'gK',
                            'EntranceAngle', 'ExitAngle'):
                    return 'Bend'
                elif hasattrs(kwargs, 'Voltage', 'Frequency', 'HarmNumber',
                              'PhaseLag', 'TimeLag'):
                    return 'RFCavity'
                elif hasattrs(kwargs, 'KickAngle'):
                    return 'Corrector'
                elif hasattrs(kwargs, 'Periodicity'):
                    return 'RingParam'
                elif hasattrs(kwargs, 'Limits'):
                    return 'Aperture'
                elif hasattrs(kwargs, 'M66'):
                    return 'M66'
                elif hasattrs(kwargs, 'GCR'):
                    return 'Monitor'
                elif hasattrs(kwargs, 'K'):
                    return 'Quadrupole'
                elif hasattrs(kwargs, 'PolynomB'):
                    try:
                        PolynomB = numpy.array(kwargs['PolynomB'],
                                               dtype=numpy.float64)
                        loworder = numpy.where(PolynomB != 0.0)[0][0]
                    except IndexError:
                        if (kwargs.get('PassMethod') ==
                            'StrMPoleSymplectic4Pass'):
                            return 'Multipole'
                        elif hasattrs(kwargs, 'BendingAngle'):
                            if float(kwargs['BendingAngle']) == 0.0:
                                return 'Drift'
                            else:
                                return 'Bend'
                        else:
                            return 'Drift'
                    if loworder == 1:
                        return 'Quadrupole'
                    elif (loworder == 2) or (float(kwargs['PolynomA'][1]) !=
                                             0.0):  # could be skew quad?
                        return 'Sextupole'
                    elif (loworder == 3) or (float(kwargs['PolynomA'][3]) !=
                                             0.0):
                        return 'Octupole'
                    else:
                        if hasattrs(kwargs, 'Length'):
                            return 'Multipole'
                        else:
                            return 'ThinMultipole'
                elif hasattrs(kwargs, 'BendingAngle'):
                    return 'Dipole'
                elif hasattrs(kwargs, 'Length'): # Don't all elems have Length?
                    if float(kwargs['Length']) != 0.0:
                        return 'Drift'
                    else:
                        return 'Marker'
                else:
                    return 'Marker'


def sanitise_class(kwargs):
    """Checks that the Class and PassMethod of the element are a valid
        combination. Some Classes and PassMethods are incompatible and would
        raise errors during calculation if left, so we raise an error here with
        a more helpful message.

    Args:
        kwargs (dict): The dictionary of keyword arguments passed to the Element
                        constructor.

    Raises:
        AttributeError: if the PassMethod and Class are incompatible.
    """
    pass_method = kwargs.get('PassMethod')
    if pass_method != None:
        pass_to_class = PASSMETHOD_MAPPING.get(pass_method)
        if (pass_method == 'IdentityPass') and (kwargs['Class'] == 'Drift'):
            kwargs['Class'] = 'Monitor'
        elif pass_to_class != None:
            if pass_to_class != kwargs['Class']:
                raise AttributeError("PassMethod {0} is not compatible with "
                                     "Class {1}.".format(pass_method,
                                                         kwargs['Class']))
        else:
            if kwargs['Class'] in ['Marker', 'Monitor', 'Drift', 'RingParam']:
                if pass_method not in ['DriftPass', 'IdentityPass']:
                    raise AttributeError("PassMethod {0} is not compatible with"
                                         " Class {1}.".format(pass_method,
                                                              kwargs['Class']))


def load_element(index, element_array):
    """
    Load what scipy produces into a pyat element object.
    """
    kwargs = {}
    kwargs['Index'] = index
    for field_name in element_array.dtype.fields:
        # Remove any surplus dimensions in arrays.
        data = numpy.squeeze(element_array[field_name])
        # Convert strings in arrays back to strings.
        if data.dtype.type is numpy.unicode_:
            data = str(data)
        kwargs[field_name] = data

    kwargs['Class'] = find_class_name(kwargs)
    sanitise_class(kwargs)

    cl = getattr(elements, kwargs['Class'])
    # Remove mandatory attributes from the keyword arguments.
    args = [kwargs.pop(attr) for attr in cl.REQUIRED_ATTRIBUTES]
    element = cl(*args, **kwargs)
    return element


def load(filename, key='RING'):
    """Load a matlab at structure into a Python at list
    """
    m = scipy.io.loadmat(filename)
    py_ring = []
    element_arrays = m[key].flat
    for i in range(len(element_arrays)):
        try:
            py_ring.append(load_element(i, element_arrays[i][0, 0]))
        except AttributeError as error_message:
            raise AttributeError('On element {0}, {1}'.format(i, error_message))
    return py_ring
