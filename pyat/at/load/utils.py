"""
Conversion utilities for creating pyat elements
"""
import os
import numpy
from warnings import warn
from distutils import sysconfig
from at import integrators
from at.lattice import elements as elt
from at.lattice.utils import AtWarning


class RingParam(elt.Element):
    """Private class for Matlab RingParam element"""
    REQUIRED_ATTRIBUTES = elt.Element.REQUIRED_ATTRIBUTES + ['Energy']
    _conversions = dict(elt.Element._conversions, Energy=float,
                        Periodicity=int)

    def __init__(self, family_name, energy, **kwargs):
        kwargs.setdefault('Periodicity', 1)
        kwargs.setdefault('PassMethod', 'IdentityPass')
        super(RingParam, self).__init__(family_name, Energy=energy, **kwargs)


# Matlab to Python class translation
_CLASS_MAP = {'drift': elt.Drift,
              'dipole': elt.Dipole, 'bend': elt.Dipole,
              'quadrupole': elt.Quadrupole, 'quad': elt.Quadrupole,
              'sextupole': elt.Sextupole, 'sext': elt.Sextupole,
              'octupole': elt.Octupole,
              'multipole': elt.Multipole,
              'thinmultipole': elt.ThinMultipole,
              'corrector': elt.Corrector,
              'rfcavity': elt.RFCavity, 'rf': elt.RFCavity,
              'monitor': elt.Monitor, 'bpm': elt.Monitor,
              'marker': elt.Marker,
              'm66': elt.M66,
              'aperture': elt.Aperture, 'ap': elt.Aperture,
              'ringparam': RingParam}

_PASS_MAP = {'DriftPass': elt.Drift,
             'BendLinearPass': elt.Dipole,
             'BndMPoleSymplectic4RadPass': elt.Dipole,
             'BndMPoleSymplectic4Pass': elt.Dipole,
             'QuadLinearPass': elt.Quadrupole,
             'CorrectorPass': elt.Corrector,
             'CavityPass': elt.RFCavity, 'RFCavityPass': elt.RFCavity,
             'ThinMPolePass': elt.ThinMultipole,
             'Matrix66Pass': elt.M66,
             'AperturePass': elt.Aperture}

# Matlab to Python attribute translation
_param_to_lattice = {'Energy': 'energy', 'Periodicity': 'periodicity',
                     'FamName': 'name'}

# Python to Matlab class translation
_matclass_map = {'Dipole': 'Bend'}

# Python to Matlab attribute translation
_matattr_map = dict(((v, k) for k, v in _param_to_lattice.items()))

# Python to Matlab type translation
_mattype_map = {int: float,
                numpy.ndarray: lambda attr: numpy.asanyarray(attr,
                                                             dtype=float)}


def hasattrs(kwargs, *attributes):
    """Check if the element would have the specified attribute(s), i.e. if they
    are in kwargs; allows checking for multiple attributes in one go.

    Args:
        kwargs (dict): The dictionary of keyword arguments passed to the
                        Element constructor.
        attributes (iterable): A list of strings, the attribute names to be
                                checked.

    Returns:
        bool: A single boolean, True if the element has any of the specified
               attributes.
    """
    for attribute in attributes:
        if attribute in kwargs:
            return True
    return False


def find_class(elem_dict, quiet=False):
    """Attempts to correctly identify the Class of the element from its kwargs.

    Args:
        elem_dict       The dictionary of keyword arguments passed to the
                        Element constructor.
    Keywords:
        quiet=False     If True, suppress the warning for non-standard classes

    Returns:
        element_class:  The guessed Class name
    """

    def low_order(key):
        polynom = numpy.array(elem_dict[key], dtype=numpy.float64).reshape(-1)
        try:
            low = numpy.where(polynom != 0.0)[0][0]
        except IndexError:
            low = -1
        return low

    class_name = elem_dict.pop('Class', '')
    try:
        return _CLASS_MAP[class_name.lower()]
    except KeyError:
        if (quiet is False) and (class_name != ''):
            warn(AtWarning("Class '{0}' does not exist.\n"
                           "{1}".format(class_name, elem_dict)))
        fam_name = elem_dict.get('FamName', '')
        try:
            return _CLASS_MAP[fam_name.lower()]
        except KeyError:
            pass_method = elem_dict.get('PassMethod', '')
            if (quiet is False) and (pass_method is ''):
                warn(AtWarning("No PassMethod provided."
                               "\n{0}".format(elem_dict)))
            elif (quiet is False) and (not pass_method.endswith('Pass')):
                warn(AtWarning("Invalid PassMethod ({0}), provided pass "
                               "methods should end in 'Pass'."
                               "\n{1}".format(pass_method, elem_dict)))
            class_from_pass = _PASS_MAP.get(pass_method)
            if class_from_pass is not None:
                return class_from_pass
            else:
                length = float(elem_dict.get('Length', 0.0))
                if hasattrs(elem_dict, 'FullGap', 'FringeInt1', 'FringeInt2',
                            'gK', 'EntranceAngle', 'ExitAngle'):
                    return elt.Dipole
                elif hasattrs(elem_dict, 'Voltage', 'Frequency', 'HarmNumber',
                              'PhaseLag', 'TimeLag'):
                    return elt.RFCavity
                elif hasattrs(elem_dict, 'Periodicity'):
                    # noinspection PyProtectedMember
                    return RingParam
                elif hasattrs(elem_dict, 'Limits'):
                    return elt.Aperture
                elif hasattrs(elem_dict, 'M66'):
                    return elt.M66
                elif hasattrs(elem_dict, 'K'):
                    return elt.Quadrupole
                elif hasattrs(elem_dict, 'PolynomB', 'PolynomA'):
                    loworder = low_order('PolynomB')
                    if loworder == 1:
                        return elt.Quadrupole
                    elif loworder == 2:
                        return elt.Sextupole
                    elif loworder == 3:
                        return elt.Octupole
                    elif (pass_method.startswith('StrMPoleSymplectic4') or
                          (length > 0)):
                        return elt.Multipole
                    else:
                        return elt.ThinMultipole
                elif hasattrs(elem_dict, 'KickAngle'):
                    return elt.Corrector
                elif length > 0.0:
                    return elt.Drift
                elif hasattrs(elem_dict, 'GCR'):
                    return elt.Monitor
                else:
                    return elt.Marker


def get_pass_method_file_name(pass_method):
    extension_list = sysconfig.get_config_vars('EXT_SUFFIX', 'SO')
    extension = set(filter(None, extension_list))
    if len(extension) == 1:
        return pass_method + extension.pop()
    else:
        return pass_method + '.so'


def element_from_dict(elem_dict, index=None, check=True, quiet=False):
    """return an AT element from a dictionary of attributes
    """

    # noinspection PyShadowingNames
    def sanitise_class(index, cls, elem_dict):
        """Checks that the Class and PassMethod of the element are a valid
            combination. Some Classes and PassMethods are incompatible and
            would raise errors during calculation if left, so we raise an error
            here with a more helpful message.

        Args:
            index:          element index
            cls:            Proposed class
            elem_dict:      he dictionary of keyword arguments passed to the
                            Element constructor.

        Raises:
            AttributeError: if the PassMethod and Class are incompatible.
        """
        def err(message, *args):
            location = ': ' if index is None else ' {0}: '.format(index)
            msg = ''.join(('Error in element', location,
                           'PassMethod {0} '.format(pass_method),
                           message.format(*args), '\n{0}'.format(elem_dict)))
            return AttributeError(msg)

        class_name = cls.__name__
        pass_method = elem_dict.get('PassMethod')
        if pass_method is not None:
            pass_to_class = _PASS_MAP.get(pass_method)
            length = float(elem_dict.get('Length', 0.0))
            file_name = get_pass_method_file_name(pass_method)
            file_path = os.path.join(integrators.__path__[0], file_name)
            if not os.path.isfile(os.path.realpath(file_path)):
                raise err("does not have a {0} file.".format(file_name))
            elif (pass_method == 'IdentityPass') and (length != 0.0):
                raise err("is not compatible with length {0}.", length)
            elif pass_to_class is not None:
                if not issubclass(cls, pass_to_class):
                    raise err("is not compatible with Class {0}.", class_name)
            elif issubclass(cls, (elt.Marker, elt.Monitor, RingParam)):
                if pass_method != 'IdentityPass':
                    raise err("is not compatible with Class {0}.", class_name)
            elif cls == elt.Drift:
                if pass_method != 'DriftPass':
                    raise err("is not compatible with Class {0}.", class_name)

    cls = find_class(elem_dict, quiet=quiet)
    if check:
        sanitise_class(index, cls, elem_dict)
    # Remove mandatory attributes from the keyword arguments.
    elem_args = (elem_dict.pop(attr, None) for attr in cls.REQUIRED_ATTRIBUTES)
    element = cls(*(arg for arg in elem_args if arg is not None), **elem_dict)
    return element


def element_to_dict(elem):
    dct = dict((k, _mattype_map.get(type(v), lambda attr: attr)(v))
               for k, v in elem.items())
    class_name = elem.__class__.__name__
    dct['Class'] = _matclass_map.get(class_name, class_name)
    return dct


def lattice_to_matlab(ring):
    dct = dict((_matattr_map.get(k, k.title()), v)
               for k, v in vars(ring).items() if not k.startswith('_'))
    famname = dct.pop('FamName')
    energy = dct.pop('Energy')
    prm = RingParam(famname, energy, **dct)
    yield element_to_dict(prm)
    for elem in ring:
        yield element_to_dict(elem)


# Kept for compatibility but should be deprecated:


CLASS_MAPPING = dict((key, cls.__name__) for (key, cls) in _CLASS_MAP.items())

PASS_MAPPING = dict((key, cls.__name__) for (key, cls) in _PASS_MAP.items())


def find_class_name(elem_dict, quiet=False):
    return find_class(elem_dict, quiet=quiet).__name__
