"""
Conversion utilities for creating pyat elements
"""
import collections
import os
import re
import numpy
from warnings import warn
from distutils import sysconfig
from at import integrators
from at.lattice import AtWarning
from at.lattice import CLASS_MAP, elements as elt
# imports necessary in' globals()' for 'eval'
# noinspection PyUnresolvedReferences
from numpy import array, uint8  # For global namespace


class RingParam(elt.Element):
    """Private class for Matlab RingParam element"""
    REQUIRED_ATTRIBUTES = elt.Element.REQUIRED_ATTRIBUTES + ['Energy',
                                                             'Periodicity']
    _conversions = dict(elt.Element._conversions, Energy=float, Periodicity=int)

    def __init__(self, family_name, energy, periodicity=1, **kwargs):
        kwargs.setdefault('Periodicity', periodicity)
        kwargs.setdefault('PassMethod', 'IdentityPass')
        super(RingParam, self).__init__(family_name, Energy=energy, **kwargs)


_alias_map = {'rbend': elt.Dipole,
              'sbend': elt.Dipole,
              'quad': elt.Quadrupole,
              'sext': elt.Sextupole,
              'rf': elt.RFCavity,
              'bpm': elt.Monitor,
              'ap': elt.Aperture,
              'ringparam': RingParam,
              'wig': elt.Wiggler}


# Matlab to Python class translation
_CLASS_MAP = dict((k.lower(), v) for k, v in CLASS_MAP.items())
_CLASS_MAP.update(_alias_map)

_PASS_MAP = {'DriftPass': elt.Drift,
             'BendLinearPass': elt.Dipole,
             'BndMPoleSymplectic4RadPass': elt.Dipole,
             'BndMPoleSymplectic4Pass': elt.Dipole,
             'QuadLinearPass': elt.Quadrupole,
             'CorrectorPass': elt.Corrector,
             'CavityPass': elt.RFCavity, 'RFCavityPass': elt.RFCavity,
             'ThinMPolePass': elt.ThinMultipole,
             'Matrix66Pass': elt.M66,
             'AperturePass': elt.Aperture,
             'GWigSymplecticPass': elt.Wiggler}

# Matlab to Python attribute translation
_param_to_lattice = {'Energy': 'energy', 'Periodicity': 'periodicity',
                     'FamName': 'name'}

# Python to Matlab class translation
_matclass_map = {'Dipole': 'Bend'}

# Python to Matlab type translation
_mattype_map = {int: float,
                numpy.ndarray: lambda attr: numpy.asanyarray(attr,
                                                             dtype=float)}

_class_to_matfunc = {
    elt.Dipole: 'atsbend',
    elt.Bend: 'atsbend',
    elt.M66: 'atM66'}


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
        if not quiet and class_name:
            warn(AtWarning("Class '{0}' does not exist.\n"
                           "{1}".format(class_name, elem_dict)))
        fam_name = elem_dict.get('FamName', '')
        try:
            return _CLASS_MAP[fam_name.lower()]
        except KeyError:
            pass_method = elem_dict.get('PassMethod', '')
            if not quiet and not pass_method:
                warn(AtWarning("No PassMethod provided."
                               "\n{0}".format(elem_dict)))
            elif not quiet and not pass_method.endswith('Pass'):
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
    ext = extension.pop() if len(extension) == 1 else '.so'
    return pass_method + ext


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
    # Create list rather than generator to ensure that elements are removed
    # from elem_dict.
    elem_args = [elem_dict.pop(attr, None) for attr in cls.REQUIRED_ATTRIBUTES]
    element = cls(*(arg for arg in elem_args if arg is not None), **elem_dict)
    return element


def element_from_string(elem_string):
    """Generate an AT-element from its repr string"""
    return eval(elem_string, globals(), CLASS_MAP)


def element_from_m(line):
    """Generate a AT-element from a line in a m-file"""

    def convert(argmt):
        """convert Matlab syntax to numpy syntax"""

        def setarray(arr):
            lns = arr.split(';')
            rr = [setarray(v) for v in lns] if len(lns) > 1 else lns[0].split()
            return '[{0}]'.format(', '.join(rr))

        if argmt.startswith('['):
            return 'array({0})'.format(setarray(argmt[1:-1]))
        else:
            return argmt

    left = line.index('(')
    right = line.rindex(')')
    cls = _CLASS_MAP[line[:left].strip()[2:]]
    arguments = [a.strip() for a in line[left + 1:right].split(',')]
    ll = len(cls.REQUIRED_ATTRIBUTES)
    if ll < len(arguments) and arguments[ll].endswith("Pass'"):
        arguments.insert(ll, "'PassMethod'")
    keys = arguments[ll::2]
    vals = arguments[ll + 1::2]
    args = [convert(v) for v in arguments[:ll]]
    keywords = ['='.join((k[1:-1], convert(v))) for k, v in zip(keys, vals)]
    elem_string = '{0}({1})'.format(cls.__name__, ', '.join(args + keywords))
    return element_from_string(elem_string)


def element_to_dict(elem):
    """Generate the Matlab structure for a AT element"""
    dct = dict((k, _mattype_map.get(type(v), lambda attr: attr)(v))
               for k, v in elem.items())
    class_name = elem.__class__.__name__
    dct['Class'] = _matclass_map.get(class_name, class_name)
    return dct


def element_to_m(elem):
    """Generate the Matlab-evaluable string for a AT element"""

    def convert(arg):
        if isinstance(arg, numpy.ndarray):
            if arg.ndim > 1:
                lns = (str(list(ln)).replace(',', '')[1:-1] for ln in arg)
                return ''.join(('[', '; '.join(lns), ']'))
            elif arg.ndim > 0:
                return str(list(arg)).replace(',', '')
            else:
                return str(arg)
        else:
            return repr(arg)

    def m_name(elclass):
        stdname = ''.join(('at', elclass.__name__.lower()))
        return _class_to_matfunc.get(elclass, stdname)

    attrs = dict(elem.items())
    args = [attrs.pop(k, getattr(elem, k)) for k in elem.REQUIRED_ATTRIBUTES]
    defelem = elem.__class__(*args)
    kwds = dict((k, v) for k, v in attrs.items()
                if not numpy.array_equal(v, getattr(defelem, k, None)))
    argstrs = [convert(arg) for arg in args]
    if 'PassMethod' in kwds:
        argstrs.append(convert(kwds.pop('PassMethod')))
    argstrs += [', '.join((repr(k), convert(v))) for k, v in kwds.items()]
    return '{0:>15}({1});...'.format(m_name(elem.__class__), ', '.join(argstrs))


# Kept for compatibility but should be deprecated:


CLASS_MAPPING = dict((key, cls.__name__) for (key, cls) in _CLASS_MAP.items())

PASS_MAPPING = dict((key, cls.__name__) for (key, cls) in _PASS_MAP.items())


def find_class_name(elem_dict, quiet=False):
    return find_class(elem_dict, quiet=quiet).__name__


def split_ignoring_parentheses(string, delimiter):
    PLACEHOLDER = "placeholder"
    substituted = string[:]
    matches = collections.deque(re.finditer("\\(.*?\\)", string))
    for match in matches:
        substituted = substituted.replace(match.group(), PLACEHOLDER, 1)
    parts = substituted.split(delimiter)
    replaced_parts = []
    for part in parts:
        if PLACEHOLDER in part:
            next_match = matches.popleft()
            part = part.replace(PLACEHOLDER, next_match.group(), 1)
        replaced_parts.append(part)
    assert not matches

    return replaced_parts
