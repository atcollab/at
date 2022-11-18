"""
Load lattices from Matlab files.
"""
from __future__ import print_function
import sys
from os.path import abspath, basename, splitext
from warnings import warn
from typing import Optional, Generator, Sequence
import scipy.io
import numpy
from ..lattice import elements, AtWarning, params_filter, AtError
from ..lattice import Element, Lattice
from .allfiles import register_format
from .utils import element_from_dict, element_from_m, RingParam
from .utils import element_to_dict, element_to_m

__all__ = ['load_mat', 'save_mat', 'load_m', 'save_m',
           'load_var']

_param_to_lattice = {'Energy': 'energy', 'Periodicity': 'periodicity',
                     'FamName': 'name', 'Particle': '_particle',
                     'cell_harmnumber': '_harmcell',
                     'HarmNumber': 'harmonic_number'}
_param_ignore = {'PassMethod', 'Length', 'cavpts'}

# Python to Matlab attribute translation
_matattr_map = dict(((v, k) for k, v in _param_to_lattice.items()))


def matfile_generator(params: dict, mat_file: str)\
        -> Generator[Element, None, None]:
    """Run through Matlab cells and generate AT elements

    Parameters:
        params:         Lattice building parameters (see :py:class:`.Lattice`)
        mat_file:       File name

    The following keys in ``params`` are used:

    ============    ===================
    **mat_key**     name of the Matlab variable containing the lattice.
                    Default: Matlab variable name if there is only one,
                    otherwise 'RING'
    **check**        Skip the coherence tests
    **quiet**       Suppress the warning for non-standard classes
    ============    ===================

    Yields:
        elem (Element): new Elements
    """
    def mclean(data):
        if data.dtype.type is numpy.str_:
            # Convert strings in arrays back to strings.
            return str(data[0]) if data.size > 0 else ''
        elif data.size == 1:
            v = data[0, 0]
            if issubclass(v.dtype.type, numpy.void):
                # Object => Return a dict
                return {f: mclean(v[f]) for f in v.dtype.fields}
            else:
                # Return a scalar
                return v
        else:
            # Remove any surplus dimensions in arrays.
            return numpy.squeeze(data)

    m = scipy.io.loadmat(params.setdefault('mat_file', mat_file))
    matvars = [varname for varname in m if not varname.startswith('__')]
    default_key = matvars[0] if (len(matvars) == 1) else 'RING'
    key = params.setdefault('mat_key', default_key)
    if key not in m.keys():
        kok = [k for k in m.keys() if '__' not in k]
        raise AtError('Selected mat_key does not exist, '
                      'please select in: {}'.format(kok))
    check = params.pop('check', True)
    quiet = params.pop('quiet', False)
    cell_array = m[key].flat
    for index, mat_elem in enumerate(cell_array):
        elem = mat_elem[0, 0]
        kwargs = {f: mclean(elem[f]) for f in elem.dtype.fields}
        yield element_from_dict(kwargs, index=index, check=check, quiet=quiet)


def ringparam_filter(params: dict, elem_iterator, *args)\
        -> Generator[Element, None, None]:
    """Run through all elements, process and optionally removes
    RingParam elements

    Parameters:
        params:         Lattice building parameters (see :py:class:`.Lattice`)
        elem_iterator:  Iterator over the lattice Elements

    Yields:
        elem (Element): new Elements

    The following keys in ``params`` are used:

    ============    ===================
    **keep_all**    keep RingParam elem_iterator as Markers
    ============    ===================

    The following keys in ``params`` are set:

    * ``name``
    * ``energy``
    * ``periodicity``
    * ``_harmnumber`` or
    * ``harmonic_number``
    * ``_particle``
    * ``_radiation``
    """
    keep_all = params.pop('keep_all', False)
    ringparams = []
    radiate = False
    for elem in elem_iterator(params, *args):
        if (elem.PassMethod.endswith('RadPass') or
                elem.PassMethod.endswith('CavityPass')):
            radiate = True
        if isinstance(elem, RingParam):
            ringparams.append(elem)
            for k, v in elem.items():
                if k not in _param_ignore:
                    params.setdefault(_param_to_lattice.get(k, k.lower()), v)
            if keep_all:
                pars = vars(elem).copy()
                name = pars.pop('FamName')
                yield elements.Marker(name, **pars)
        else:
            yield elem
    params['_radiation'] = radiate

    if len(ringparams) > 1:
        warn(AtWarning('More than 1 RingParam element, the 1st one is used'))


def load_mat(filename: str, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from a Matlab mat-file

    Parameters:
        filename:           Name of a '.mat' file

    Keyword Args:
        mat_key (str):      Name of the Matlab variable containing
          the lattice. Default: Matlab variable name if there is only one,
          otherwise 'RING'
        check (bool):       Run the coherence tests. Default:
          :py:obj:`True`
        quiet (bool):       Suppress the warning for non-standard
          classes. Default: :py:obj:`False`
        keep_all (bool):    Keep RingParam elements as Markers.
          Default: :py:obj:`False`
        name (str):         Name of the lattice. Default: taken from
          the lattice, or ``''``
        energy (float):     Energy of the lattice [eV]. Default: taken
          from the lattice elements
        periodicity(int):   Number of periods. Default: taken from the
          elements, or 1
        *:                  All other keywords will be set as Lattice
          attributes

    Returns:
        lattice (Lattice):  New :py:class:`.Lattice` object

    See Also:
        :py:func:`.load_lattice` for a generic lattice-loading function.
   """
    if 'key' in kwargs:  # process the deprecated 'key' keyword
        kwargs.setdefault('mat_key', kwargs.pop('key'))
    return Lattice(ringparam_filter, matfile_generator, abspath(filename),
                   iterator=params_filter, **kwargs)


def mfile_generator(params: dict, m_file: str)\
        -> Generator[Element, None, None]:
    """Run through the lines of a Matlab m-file and generate AT elements

    Parameters:
        params:         Lattice building parameters (see :py:class:`.Lattice`)
        m_file:         File name

    Yields:
        elem (Element): new Elements
"""
    with open(params.setdefault('m_file', m_file), 'rt') as file:
        _ = next(file)  # Matlab function definition
        _ = next(file)  # Cell array opening
        for lineno, line in enumerate(file):
            if line.startswith('};'):
                break
            try:
                elem = element_from_m(line)
            except ValueError:
                warn(AtWarning('Invalid line {0} skipped.'.format(lineno)))
                continue
            except KeyError:
                warn(AtWarning('Line {0}: Unknown class.'.format(lineno)))
                continue
            else:
                yield elem


def load_m(filename: str, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from a Matlab m-file

    Parameters:
        filename:           Name of a '.m' file

    Keyword Args:
        keep_all (bool):    Keep RingParam elements as Markers.
          Default: :py:obj:`False`
        name (str):         Name of the lattice. Default: taken from
          the lattice, or ``''``
        energy (float):     Energy of the lattice [eV]. Default: taken
          from the lattice elements
        periodicity(int):   Number of periods. Default: taken from the
          elements, or 1
        *:                  All other keywords will be set as Lattice
          attributes

    Returns:
        lattice (Lattice):  New :py:class:`.Lattice` object

    See Also:
        :py:func:`.load_lattice` for a generic lattice-loading function.
    """
    return Lattice(ringparam_filter, mfile_generator, abspath(filename),
                   iterator=params_filter, **kwargs)


def load_var(matlat: Sequence[dict], **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice` from a Matlab cell array

    Parameters:
        matlat:             Matlab lattice

    Keyword Args:
        keep_all (bool):    Keep RingParam elements as Markers.
          Default: :py:obj:`False`
        name (str):         Name of the lattice. Default: taken from
          the lattice, or ``''``
        energy (float):     Energy of the lattice [eV]. Default: taken
          from the lattice elements
        periodicity(int):   Number of periods. Default: taken from the
          elements, or 1
        *:                  All other keywords will be set as Lattice
          attributes

    Returns:
        lattice (Lattice):  New :py:class:`.Lattice` object
    """
    # noinspection PyUnusedLocal
    def var_generator(params, latt):
        for elem in latt:
            yield element_from_dict(elem)

    return Lattice(ringparam_filter, var_generator, matlat,
                   iterator=params_filter, **kwargs)


def matlab_ring(ring) -> Generator[Element, None, None]:
    """Prepend a RingParam element to a lattice"""
    dct = dict((_matattr_map.get(k, k.title()), v)
               for k, v in ring.attrs.items())
    famname = dct.pop('FamName')
    energy = dct.pop('Energy')
    yield RingParam(famname, energy, **dct)
    for elem in ring:
        yield elem


def save_mat(ring: Lattice, filename: str,
             mat_key: str = 'RING') -> None:
    """Save a :py:class:`.Lattice` as a Matlab mat-file

    Parameters:
        ring:           Lattice description
        filename:       Name of the '.mat' file
        mat_key (str):  Name of the Matlab variable containing
          the lattice. Default: ``'RING'``

    See Also:
        :py:func:`.save_lattice` for a generic lattice-saving function.
    """
    lring = tuple((element_to_dict(elem),) for elem in matlab_ring(ring))
    scipy.io.savemat(filename, {mat_key: lring}, long_field_names=True)


def save_m(ring: Lattice, filename: Optional[str] = None) -> None:
    """Save a :py:class:`.Lattice` as a Matlab m-file

    Parameters:
        ring:           Lattice description
        filename:       Name of the '.m' file. Default: outputs on
          :py:obj:`sys.stdout`

    See Also:
        :py:func:`.save_lattice` for a generic lattice-saving function.
    """

    def save(file):
        print('ring = {...', file=file)
        for elem in matlab_ring(ring):
            print(element_to_m(elem), file=file)
        print('};', file=file)

    if filename is None:
        save(sys.stdout)
    else:
        with open(filename, 'wt') as mfile:
            [funcname, _] = splitext(basename(filename))
            print('function ring = {0}()'.format(funcname), file=mfile)
            save(mfile)
            print('end', file=mfile)


register_format('.mat', load_mat, save_mat, descr='Matlab binary mat-file')
register_format('.m', load_m, save_m, descr='Matlab text m-file')
