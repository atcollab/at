"""
Load lattices from Matlab files.
"""
from __future__ import print_function
import sys
from os.path import abspath
from warnings import warn
import scipy.io
import numpy
from at.lattice import elements, Lattice, AtWarning, params_filter
from at.load import register_format
from at.load.utils import element_from_dict, element_from_m, RingParam
from at.load.utils import element_to_dict, element_to_m

__all__ = ['ringparam_filter', 'load_mat', 'save_mat', 'load_m', 'save_m',
           'load_var']

_param_to_lattice = {'Energy': 'energy', 'Periodicity': 'periodicity',
                     'FamName': 'name'}
_param_ignore = {'PassMethod', 'Length'}

# Python to Matlab attribute translation
_matattr_map = dict(((v, k) for k, v in _param_to_lattice.items()))


def matfile_generator(params, mat_file):
    """
    run through matlab cells and generate AT elements

    KEYWORDS
        mat_file        name of the .mat file
        mat_key         name of the Matlab variable containing the lattice.
                        Default: Matlab variable name if there is only one,
                        otherwise 'RING'
        check=True      if False, skip the coherence tests
        quiet=False     If True, suppress the warning for non-standard classes
    """
    m = scipy.io.loadmat(params.setdefault('mat_file', mat_file))
    matvars = [varname for varname in m if not varname.startswith('__')]
    default_key = matvars[0] if (len(matvars) == 1) else 'RING'
    key = params.setdefault('mat_key', default_key)
    check = params.pop('check', True)
    quiet = params.pop('quiet', False)
    cell_array = m[key].flat
    for index, mat_elem in enumerate(cell_array):
        element_array = mat_elem[0, 0]
        kwargs = {}
        for field_name in element_array.dtype.fields:
            # Remove any surplus dimensions in arrays.
            data = numpy.squeeze(element_array[field_name])
            # Convert strings in arrays back to strings.
            if data.dtype.type is numpy.unicode_:
                data = str(data)
            kwargs[field_name] = data

        yield element_from_dict(kwargs, index=index, check=check, quiet=quiet)


def ringparam_filter(params, elem_iterator, *args):
    """"
    Run through all elements, process and optionally removes RingParam elements

    KEYWORDS
        keep_all=False  if True, keep RingParam elem_iterator as Markers
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


def load_mat(filename, **kwargs):
    """Create a lattice object from a matmab mat-file

    PARAMETERS
        filename        name of a '.mat' file

    KEYWORDS
        mat_key         name of the Matlab variable containing the lattice.
                        Default: Matlab variable name if there is only one,
                        otherwise 'RING'
        check=True      if False, skip the coherence tests
        quiet=False     If True, suppress the warning for non-standard classes
        keep_all=False  if True, keep RingParam elements as Markers
        name            Name of the lattice
                        (default: taken from the lattice, or '')
        energy          Energy of the lattice
                        (default: taken from the elements)
        periodicity     Number of periods
                        (default: taken from the elements, or 1)
        *               all other keywords will be set as Lattice attributes

    OUTPUT
        Lattice object
    """
    if 'key' in kwargs:  # process the deprecated 'key' keyword
        kwargs.setdefault('mat_key', kwargs.pop('key'))
    return Lattice(ringparam_filter, matfile_generator, abspath(filename),
                   iterator=params_filter, **kwargs)


def mfile_generator(params, m_file):
    """Run through all lines of a m-file and generates AT elements"""
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
                warn(AtWarning('Line {0}: Unknown class'))
                continue
            else:
                yield elem


def load_m(filename, **kwargs):
    """Create a lattice object from a matlab m-file

    PARAMETERS
        filename        name of a '.m' file

    KEYWORDS
        keep_all=False  if True, keep RingParam elements as Markers
        name            Name of the lattice
                        (default: taken from the elements, or '')
        energy          Energy of the lattice
                        (default: taken from the elements)
        periodicity     Number of periods
                        (default: taken from the elements, or 1)
        *               all other keywords will be set as Lattice attributes

    OUTPUT
        Lattice object
    """
    return Lattice(ringparam_filter, mfile_generator, abspath(filename),
                   iterator=params_filter, **kwargs)


def load_var(matlat, **kwargs):
    """Create a lattice from a Matlab cell array"""

    # noinspection PyUnusedLocal
    def var_generator(params, latt):
        for elem in latt:
            yield element_from_dict(elem)

    return Lattice(ringparam_filter, var_generator, matlat,
                   iterator=params_filter, **kwargs)


def matlab_ring(ring):
    """Prepend a RingParam element to a lattice"""
    dct = dict((_matattr_map.get(k, k.title()), v)
               for k, v in vars(ring).items() if not k.startswith('_'))
    famname = dct.pop('FamName')
    energy = dct.pop('Energy')
    yield RingParam(famname, energy, **dct)
    for elem in ring:
        yield elem


def save_mat(ring, filename, mat_key='RING'):
    """Save a Lattice object as a Matlab mat-file

    PARAMETERS
        ring            Lattice object
        filename        name of the '.mat' file

    KEYWORDS
        mat_key='RING'  Name of the Matlab variable representing the lattice
    """
    lring = tuple((element_to_dict(elem),) for elem in matlab_ring(ring))
    scipy.io.savemat(filename, {mat_key: lring})


def save_m(ring, filename=None):
    """Save a lattice as a Matlab m-file

    PARAMETERS
        ring            Lattice object
        filename=None   name of the '.m' file. Default: output on sys.stdout
    """

    def save(file):
        print('function ring = {0}()'.format(ring.name), file=file)
        print('ring = {...', file=file)
        for elem in matlab_ring(ring):
            print(element_to_m(elem), file=file)
        print('};', file=file)
        print('end', file=file)

    if filename is None:
        save(sys.stdout)
    else:
        with open(filename, 'wt') as mfile:
            save(mfile)


register_format('.mat', load_mat, save_mat, descr='Matlab binary mat-file')
register_format('.m', load_m, save_m, descr='Matlab text m-file')
