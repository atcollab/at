"""
Conversion utilities for creating pyat elements
"""
import collections
import os
import re
import numpy
from warnings import warn
from typing import Optional
import sysconfig
from at import integrators
from at.lattice import AtWarning
from at.lattice import CLASS_MAP, elements as elt
from at.lattice import idtable_element
from at.lattice import Particle, Element
# imports necessary in' globals()' for 'eval'
# noinspection PyUnresolvedReferences
from numpy import array, uint8  # For global namespace

_ext_suffix = sysconfig.get_config_var('EXT_SUFFIX')


def _particle(value):
    if isinstance(value, Particle):
        # Create from python: save_mat
        return value
    else:
        # Create from Matlab: load_mat
        name = value.pop('name')
        return Particle(name, **value)


class RingParam(elt.Element):
    """Private class for Matlab RingParam element

    :meta private:
    """
    # noinspection PyProtectedMember
    _BUILD_ATTRIBUTES = elt.Element._BUILD_ATTRIBUTES + ['Energy',
                                                         'Periodicity']
    _conversions = dict(elt.Element._conversions, Energy=float,
                        Periodicity=int, Particle=_particle)

    def __init__(self, family_name, energy, periodicity=1, **kwargs):
        kwargs.setdefault('Energy', energy)
        kwargs.setdefault('Periodicity', periodicity)
        kwargs.setdefault('PassMethod', 'IdentityPass')
        super(RingParam, self).__init__(family_name, **kwargs)


_alias_map = {'rbend': elt.Dipole,
              'sbend': elt.Dipole,
              'quad': elt.Quadrupole,
              'sext': elt.Sextupole,
              'rf': elt.RFCavity,
              'bpm': elt.Monitor,
              'ap': elt.Aperture,
              'ringparam': RingParam,
              'wig': elt.Wiggler,
              'insertiondevicekickmap': idtable_element.InsertionDeviceKickMap
              }


# Matlab to Python class translation
_CLASS_MAP = dict((k.lower(), v) for k, v in CLASS_MAP.items())
_CLASS_MAP.update(_alias_map)

_PASS_MAP = {'BendLinearPass': elt.Dipole,
             'BndMPoleSymplectic4RadPass': elt.Dipole,
             'BndMPoleSymplectic4Pass': elt.Dipole,
             'QuadLinearPass': elt.Quadrupole,
             'StrMPoleSymplectic4Pass': elt.Multipole,
             'StrMPoleSymplectic4RadPass': elt.Multipole,
             'CorrectorPass': elt.Corrector,
             'CavityPass': elt.RFCavity, 'RFCavityPass': elt.RFCavity,
             'ThinMPolePass': elt.ThinMultipole,
             'Matrix66Pass': elt.M66,
             'AperturePass': elt.Aperture,
             'IdTablePass': idtable_element.InsertionDeviceKickMap,
             'GWigSymplecticPass': elt.Wiggler}

# Matlab to Python attribute translation
_param_to_lattice = {'Energy': 'energy', 'Periodicity': 'periodicity',
                     'FamName': 'name'}

# Python to Matlab class translation
_matclass_map = {
        'Dipole': 'Bend',
        'InsertionDeviceKickMap': 'InsertionDeviceKickMap'
        }

# Python to Matlab type translation
_mattype_map = {int: float,
                numpy.ndarray: lambda attr: numpy.asanyarray(attr),
                Particle: lambda attr: attr.to_dict()}

_class_to_matfunc = {
    elt.Dipole: 'atsbend',
    elt.Bend: 'atsbend',
    elt.M66: 'atM66',
    idtable_element.InsertionDeviceKickMap: 'atinsertiondevicekickmap'
    }


def hasattrs(kwargs: dict, *attributes) -> bool:
    """Checks the presence of keys in a :py:class:`dict`

    Returns :py:obj:`True` if any of the ``attributes`` is in ``kwargs``

    Args:
        kwargs:     The dictionary of keyword arguments passed to the
          Element constructor.
        attributes: A list of strings, the attribute names to be checked.

    Returns:
        found (bool):   :py:obj:`True` if the element has any of the specified
          attributes.
    """
    for attribute in attributes:
        if attribute in kwargs:
            return True
    return False


def find_class(elem_dict: dict, quiet: bool = False) -> type(Element):
    """Identify the class of an element from its attributes

    Args:
        elem_dict:      The dictionary of keyword arguments passed to the
                        Element constructor.
        quiet:          Suppress the warning for non-standard classes

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
            class_doesnotexist_warning = ("Class '{0}' does not exist.\n"
                                          "{1}".format(class_name, elem_dict))
            warn(AtWarning(class_doesnotexist_warning))
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
                elif pass_method == 'IdentityPass':
                    return elt.Marker
                else:
                    return elt.Element


def element_from_dict(elem_dict: dict, index: Optional[int] = None,
                      check: bool = True, quiet: bool = False) -> Element:
    """Builds an :py:class:`.Element` from a dictionary of attributes

    Parameters:
        elem_dict:      Dictionary of element attributes
        index:          Element index
        check:          Check the compatibility of class and PassMethod
        quiet:          Suppress the warning for non-standard classes

    Returns:
        elem (Element): new :py:class:`.Element`
    """

    # noinspection PyShadowingNames
    def sanitise_class(index, cls, elem_dict):
        """Checks that the Class and PassMethod of the element are a valid
            combination. Some Classes and PassMethods are incompatible and
            would raise errors during calculation, so we send a
            warning here.

        Args:
            index:          element index
            cls:            Proposed class
            elem_dict:      The dictionary of keyword arguments passed to the
                            Element constructor.

        Raises:
            AttributeError: if the PassMethod and Class are incompatible.
        """
        def err(message, *args):
            location = ': ' if index is None else ' {0}: '.format(index)
            msg = ''.join(('Error in element', location,
                           'PassMethod {0} '.format(pass_method),
                           message.format(*args), '\n{0}'.format(elem_dict)))
            return AtWarning(msg)

        class_name = cls.__name__
        pass_method = elem_dict.get('PassMethod')
        if pass_method is not None:
            pass_to_class = _PASS_MAP.get(pass_method)
            length = float(elem_dict.get('Length', 0.0))
            file_name = pass_method + _ext_suffix
            file_path = os.path.join(integrators.__path__[0], file_name)
            if not os.path.isfile(os.path.realpath(file_path)):
                warn(err(" is missing {0}.".format(file_name)))
            elif (pass_method == 'IdentityPass') and (length != 0.0):
                warn(err("is not compatible with length {0}.", length))
            elif pass_to_class is not None:
                if not issubclass(cls, pass_to_class):
                    warn(err("is not compatible with Class {0}.", class_name))

    cls = find_class(elem_dict, quiet=quiet)
    if check:
        sanitise_class(index, cls, elem_dict)
    # Remove mandatory attributes from the keyword arguments.
    # Create list rather than generator to ensure that elements are removed
    # from elem_dict.
    elem_args = [elem_dict.pop(attr, None) for attr in cls._BUILD_ATTRIBUTES]
    element = cls(*(arg for arg in elem_args if arg is not None), **elem_dict)
    return element


def element_from_string(elem_string: str) -> Element:
    """Builds an :py:class:`.Element` from its python :py:func:`repr` string

    Parameters:
        elem_string:    String representation of an :py:class:`.Element`

    Returns:
        elem (Element): new :py:class:`.Element`
    """
    return eval(elem_string, globals(), CLASS_MAP)


def element_from_m(line: str) -> Element:
    """Builds an :py:class:`.Element` from a line in an m-file

    Parameters:
        line:           Matlab string representation of an :py:class:`.Element`

    Returns:
        elem (Element): new :py:class:`.Element`
    """
    def argsplit(value):
        return [a.strip() for a in split_ignoring_parentheses(value, ',')]

    def makedir(mat_struct):
        """Build directory from Matlab struct arguments"""
        def pairs(it):
            while True:
                try:
                    a = next(it)
                except StopIteration:
                    break
                yield eval(a), convert(next(it))
        return dict(pairs(iter(mat_struct)))

    def makearray(mat_arr):
        """Build numpy array for Matlab array syntax"""
        def arraystr(arr):
            lns = arr.split(';')
            rr = [arraystr(v) for v in lns] if len(lns) > 1 else lns[0].split()
            return '[{0}]'.format(', '.join(rr))
        return eval('numpy.array({0})'.format(arraystr(mat_arr)))

    def convert(value):
        """convert Matlab syntax to numpy syntax"""
        if value.startswith('['):
            result = makearray(value[1:-1])
        elif value.startswith('struct'):
            result = makedir(argsplit(value[7:-1]))
        else:
            result = eval(value)
        return result

    left = line.index('(')
    right = line.rindex(')')
    matcls = line[:left].strip()[2:]
    cls = _CLASS_MAP[matcls]
    arguments = argsplit(line[left + 1:right])
    ll = len(cls._BUILD_ATTRIBUTES)
    if ll < len(arguments) and arguments[ll].endswith("Pass'"):
        arguments.insert(ll, "'PassMethod'")
    args = [convert(v) for v in arguments[:ll]]
    kwargs = makedir(arguments[ll:])
    if matcls == 'rbend':
        # the Matlab 'rbend' has no equivalent in PyAT. This adds parameters
        # necessary for using the python sector bend
        halfangle = 0.5 * args[2]
        kwargs.setdefault('EntranceAngle', halfangle)
        kwargs.setdefault('ExitAngle', halfangle)
    return cls(*args, **kwargs)


def element_to_dict(elem: Element) -> dict:
    """Builds the Matlab structure of an :py:class:`.Element`

    Parameters:
        elem:           :py:class:`.Element`

    Returns:
        dct (dict):     Dictionary of :py:class:`.Element` attributes
    """
    dct = dict((k, _mattype_map.get(type(v), lambda attr: attr)(v))
               for k, v in elem.items())
    class_name = elem.__class__.__name__
    dct['Class'] = _matclass_map.get(class_name, class_name)
    return dct


def element_to_m(elem: Element) -> str:
    """Builds the Matlab-evaluable string for an :py:class:`.Element`

    Parameters:
        elem:           :py:class:`.Element`

    Returns:
        mstr (str):     Matlab string representation of the
          :py:class:`.Element` attributes
    """

    def convert(arg):
        def convert_dict(pdir):
            def scan(d):
                for k, v in d.items():
                    yield convert(k)
                    yield convert(v)
            return 'struct({0})'.format(', '.join(scan(pdir)))

        def convert_array(arr):
            if arr.ndim > 1:
                lns = (str(list(ln)).replace(',', '')[1:-1] for ln in arr)
                return ''.join(('[', '; '.join(lns), ']'))
            elif arr.ndim > 0:
                return str(list(arr)).replace(',', '')
            else:
                return str(arr)

        if isinstance(arg, numpy.ndarray):
            return convert_array(arg)
        elif isinstance(arg, dict):
            return convert_dict(arg)
        elif isinstance(arg, Particle):
            return convert_dict(arg.to_dict())
        else:
            return repr(arg)

    def m_name(elclass):
        stdname = ''.join(('at', elclass.__name__.lower()))
        return _class_to_matfunc.get(elclass, stdname)

    attrs = dict(elem.items())
    # noinspection PyProtectedMember
    args = [attrs.pop(k, getattr(elem, k)) for k in elem._BUILD_ATTRIBUTES]
    defelem = elem.__class__(*args)
    kwds = dict((k, v) for k, v in attrs.items()
                if not numpy.array_equal(v, getattr(defelem, k, None)))
    argstrs = [convert(arg) for arg in args]
    if 'PassMethod' in kwds:
        argstrs.append(convert(kwds.pop('PassMethod')))
    argstrs += [', '.join((repr(k), convert(v))) for k, v in kwds.items()]
    return '{0:>15}({1});...'.format(m_name(elem.__class__),
                                     ', '.join(argstrs))


# Kept for compatibility but should be deprecated:


CLASS_MAPPING = dict((key, cls.__name__) for (key, cls) in _CLASS_MAP.items())

PASS_MAPPING = dict((key, cls.__name__) for (key, cls) in _PASS_MAP.items())


def find_class_name(elem_dict, quiet=False):
    """Derive the class name of an Element from its attributes"""
    return find_class(elem_dict, quiet=quiet).__name__


def split_ignoring_parentheses(string, delimiter):
    placeholder = "placeholder"
    substituted = string[:]
    matches = collections.deque(re.finditer("\\(.*?\\)", string))
    for match in matches:
        substituted = substituted.replace(match.group(), placeholder, 1)
    parts = substituted.split(delimiter)
    replaced_parts = []
    for part in parts:
        if placeholder in part:
            next_match = matches.popleft()
            part = part.replace(placeholder, next_match.group(), 1)
        replaced_parts.append(part)
    assert not matches

    return replaced_parts
