"""
Load lattices from Matlab files.
"""

from __future__ import annotations

__all__ = ["load_mat", "save_mat", "load_m", "save_m", "load_var"]

import os
import sys
from collections.abc import Generator, Sequence
from math import isfinite
from os.path import abspath, basename, splitext
from typing import Any
from warnings import warn

import h5py
import numpy as np
import scipy.io

# imports necessary in 'globals()' for 'eval'
from numpy import array  # noqa: F401
from numpy import nan as NaN  # noqa: F401
from numpy import uint8  # noqa: F401

from ..lattice import (
    AtError,
    AtWarning,
    Element,
    Filter,
    Lattice,
    Particle,
    elements,
    params_filter,
)
from .allfiles import register_format
from .utils import (
    _CLASS_MAP,
    RingParam,
    _drop_attrs,
    keep_elements,
    split_ignoring_parentheses,
)

# Translation of RingParam attributes
_m2p = {
    "FamName": "name",
    "Energy": "energy",
    "Periodicity": "periodicity",
    "Particle": "particle",
    "cell_harmnumber": "cell_harmnumber",  # necessary: property
    "beam_current": "beam_current",  # necessary: property
    "PassMethod": None,  # Useless Matlab attributes
    "Length": None,
    "cavpts": None,
    "Mat_File": None,  # These are erroneous attributes saved in
    "Mat_Key": None,  # RingParam by old versions
    "Beam_Current": None,
    "Nbunch": None,
}
_p2m = {v: k for k, v in _m2p.items() if v is not None}
# Attribute to drop when writing a file
_p2m.update(_drop_attrs)

# Matlab constructor function
# Default: "".join(("at", element_class.__name__.lower()))
_mat_constructor = {
    "Dipole": "atsbend",
    "M66": "atM66",
}


def _mat_encoder(v):
    """type encoding for .mat files"""
    if isinstance(v, int):
        return float(v)
    if isinstance(v, Particle):
        return v.to_dict()
    if isinstance(v, (list, tuple)):
        return np.array(v, dtype=object)
    return v


def _matfile_generator(
    params: dict[str, Any], mat_file: str
) -> Generator[Element, None, None]:
    """Run through Matlab cells and generate AT elements.

    Arguments:
        params: parameter dictionary
        mat_file: matlab file name

    Yields:
        pyat Element from dictionary
    """

    def mclean(data: any) -> any:
        if data.dtype.type is np.str_:
            # Convert strings in arrays back to strings.
            dataout = str(data[0]) if data.size > 0 else ""
        elif data.size == 1:
            vdata = data[0, 0]
            if issubclass(vdata.dtype.type, np.void):
                # Object => Return a dict
                dataout = {f: mclean(vdata[f]) for f in vdata.dtype.fields}
            else:
                # Return a scalar
                dataout = vdata
        else:
            # Remove any surplus dimensions in arrays.
            dataout = np.squeeze(data)
        return dataout

    def mcleanhdf5(data: any) -> any:
        matlab_class = data.attrs["MATLAB_class"]
        if matlab_class == b"struct":
            # Return a dict from recursion
            dataout = {f: mcleanhdf5(data[f]) for f in data.keys()}
        elif matlab_class == b"char":
            # Convert to string
            dataout = "".join(chr(i) for i in np.asarray(data).flatten())
        else:
            # e.g. matlab_class == b"double":
            # Remove any surplus dimensions in arrays.
            dataout = np.squeeze(np.asarray(data))
        return dataout

    def define_default_key(
        params: dict, mat_input: any, ignore_chars: str = ""
    ) -> tuple:
        matvars = [
            varname for varname in mat_input if not varname.startswith(ignore_chars)
        ]
        default_key = matvars[0] if (len(matvars) == 1) else "RING"
        key = params.setdefault("use", default_key)
        if key not in mat_input.keys():
            kok = [k for k in mat_input.keys() if "__" not in k]
            raise AtError(
                f"Selected '{key}' variable does not exist, please select in: {kok}"
            )
        return params, key

    # noinspection PyUnresolvedReferences
    check = params.pop("check", True)
    quiet = params.pop("quiet", False)
    matlabfile_ver = scipy.io.matlab.matfile_version(mat_file)
    if matlabfile_ver < (2, 0):
        mat_input = scipy.io.loadmat(params.setdefault("in_file", mat_file))
        params, key = define_default_key(params, mat_input, ignore_chars="__")
        cell_array = mat_input[key].flat
        for index, mat_elem in enumerate(cell_array):
            elem = mat_elem[0, 0]
            kwargs = {f: mclean(elem[f]) for f in elem.dtype.fields}
            yield Element.from_matlab(kwargs, index=index, check=check, quiet=quiet)
    else:
        mat_input = h5py.File(mat_file)
        params, key = define_default_key(params, mat_input, ignore_chars="#")
        cell_array = mat_input[key][0]
        for index, ref_elem in enumerate(cell_array):
            elem = mat_input[ref_elem]
            kwargs = {f: mcleanhdf5(elem[f]) for f in elem.keys()}
            yield Element.from_matlab(kwargs, index=index, check=check, quiet=quiet)


def ringparam_filter(
    params: dict[str, Any], elem_iterator: Filter, *args
) -> Generator[Element, None, None]:
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
    keep_all = params.pop("keep_all", False)
    ringparams = []
    radiate = False
    for elem in elem_iterator(params, *args):
        if elem.PassMethod.endswith("RadPass") or elem.PassMethod.endswith(
            "CavityPass"
        ):
            radiate = True
        if isinstance(elem, RingParam):
            ringparams.append(elem)
            for k, v in elem.items():
                k2 = _m2p.get(k, k)
                if k2 is not None:
                    if k2 != "cell_harmnumber" or isfinite(v):
                        params.setdefault(k2, v)
            if keep_all:
                pars = vars(elem).copy()
                name = pars.pop("FamName")
                yield elements.Marker(name, tag="RingParam", **pars)
        else:
            yield elem
    params["_radiation"] = radiate
    params.setdefault("periodicity", 0)

    if len(ringparams) > 1:
        warn(
            AtWarning("More than 1 RingParam element, the 1st one is used"),
            stacklevel=2,
        )


def load_mat(filename: str, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from a Matlab mat-file

    Parameters:
        filename:           Name of a '.mat' file

    Keyword Args:
        use (str):          Name of the Matlab variable containing
          the lattice. Default: it there is a single variable, use it, otherwise
          select 'RING'
        mat_key (str):      deprecated alias for *use*
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
    if "key" in kwargs:  # process the deprecated 'key' keyword
        kwargs.setdefault("use", kwargs.pop("key"))
    if "mat_key" in kwargs:  # process the deprecated 'mat_key' keyword
        kwargs.setdefault("use", kwargs.pop("mat_key"))
    return Lattice(
        ringparam_filter,
        _matfile_generator,
        abspath(filename),
        iterator=params_filter,
        **kwargs,
    )


def _element_from_m(line: str) -> Element:
    """Builds an :py:class:`.Element` from a line in an m-file

    Parameters:
        line:           Matlab string representation of an :py:class:`.Element`

    Returns:
        elem (Element): new :py:class:`.Element`
    """

    def argsplit(value):
        return [a.strip() for a in split_ignoring_parentheses(value)]

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

    def makearray(mat_arr: str) -> np.ndarray:
        """Build numpy array from Matlab array as strings.

        Accepts matlab 1D or 2D array repr [...], [...;...], or [[...];[...]].

        Arguments:
            mat_arr: matlab style array.

        Returns:
            numpy style ndarray.
        """
        lns = mat_arr.replace("[", "").replace("]", "").split(";")
        return np.squeeze(np.array([row.split() for row in lns], dtype=np.float64))

    def convert(value):
        """convert Matlab syntax to numpy syntax"""
        if value.startswith("["):
            result = makearray(value)
        elif value.startswith("struct"):
            result = makedir(argsplit(value[7:-1]))
        else:
            result = eval(value)
        return result

    left = line.index("(")
    right = line.rindex(")")
    matcls = line[:left].strip()[2:]
    cls = _CLASS_MAP[matcls]
    arguments = argsplit(line[left + 1 : right])
    ll = len(cls._BUILD_ATTRIBUTES)
    if ll < len(arguments) and arguments[ll].endswith("Pass'"):
        arguments.insert(ll, "'PassMethod'")
    args = [convert(v) for v in arguments[:ll]]
    kwargs = makedir(arguments[ll:])
    if matcls == "rbend":
        # the Matlab 'rbend' has no equivalent in PyAT. This adds parameters
        # necessary for using the python sector bend
        halfangle = 0.5 * args[2]
        kwargs.setdefault("EntranceAngle", halfangle)
        kwargs.setdefault("ExitAngle", halfangle)
    return cls(*args, **kwargs)


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

    def mfile_generator(params: dict, m_file: str) -> Generator[Element, None, None]:
        """Run through the lines of a Matlab m-file and generate AT elements"""
        with open(params.setdefault("in_file", m_file)) as file:
            _ = next(file)  # Matlab function definition
            _ = next(file)  # Cell array opening
            for lineno, line in enumerate(file):
                if line.startswith("};"):
                    break
                try:
                    elem = _element_from_m(line)
                except ValueError:
                    warn(AtWarning(f"Invalid line {lineno} skipped."), stacklevel=2)
                    continue
                except KeyError:
                    warn(AtWarning(f"Line {lineno}: Unknown class."), stacklevel=2)
                    continue
                else:
                    yield elem

    return Lattice(
        ringparam_filter,
        mfile_generator,
        abspath(filename),
        iterator=params_filter,
        **kwargs,
    )


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
            yield Element.from_matlab(elem)

    return Lattice(
        ringparam_filter, var_generator, matlat, iterator=params_filter, **kwargs
    )


def matlab_ring(ring: Lattice) -> Generator[Element, None, None]:
    """Prepend a RingParam element to a lattice"""

    def required(rng):
        # Public lattice attributes
        params = {k: v for k, v in vars(rng).items() if not k.startswith("_")}
        # Output the required attributes/properties
        for kp, km in _p2m.items():
            try:
                v = getattr(rng, kp)
            except AttributeError:
                pass
            else:
                params.pop(kp, None)
                if km is not None:
                    yield km, v
        # Output the remaining attributes
        yield from params.items()

    dct = dict(required(ring))
    yield RingParam(**dct)
    yield from keep_elements(ring)


def save_mat(ring: Lattice, filename: str, **kwargs) -> None:
    """Save a :py:class:`.Lattice` as a Matlab mat-file

    Parameters:
        ring:       Lattice description
        filename:   Name of the '.mat' file

    Keyword Args:
        use (str):  Name of the Matlab variable containing the lattice, Default: "RING"
        mat_key (str): Deprecated, alias for *use*

    See Also:
        :py:func:`.save_lattice` for a generic lattice-saving function.
    """
    # Ensure the lattice is a Matlab column vector: list(list)
    use = kwargs.pop("mat_key", "RING")  # For backward compatibility
    use = kwargs.pop("use", use)
    lring = [el.to_matlab(encoder=_mat_encoder) for el in matlab_ring(ring)]
    scipy.io.savemat(filename, {use: lring}, long_field_names=True)


def _element_to_m(elem: Element) -> str:
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

            return "struct({})".format(", ".join(scan(pdir)))

        def convert_array(arr):
            max_array = max(1000, np.prod(arr.shape)) + 1
            mod_opt = {"threshold": max_array, "max_line_width": np.inf}
            if arr.ndim > 1:
                # replace endline character by ; to indicate the end of a 1D array
                # replace [SPACE by [ to remove extra space from +sign of
                # first array element
                return (
                    np.array2string(arg, **mod_opt)
                    .replace("\n", ";")
                    .replace("[ ", "[")
                )
            elif arr.ndim > 0:
                return np.array2string(arg, **mod_opt).replace("[ ", "[")
            else:
                return str(arr)

        def convert_list(lst):
            return f"{{{{{str(lst)[1:-1]}}}}}"

        if isinstance(arg, np.ndarray):
            return convert_array(arg)
        elif isinstance(arg, np.number):
            return str(arg)
        elif isinstance(arg, dict):
            return convert_dict(arg)
        elif isinstance(arg, (list, tuple)):
            return convert_list(arg)
        elif isinstance(arg, Particle):
            return convert_dict(arg.to_dict())
        else:
            return repr(arg)

    def m_name(elclass):
        classname = elclass.__name__
        return _mat_constructor.get(classname, "".join(("at", classname.lower())))

    _, args, kwds = elem.definition
    argstrs = [convert(arg) for arg in args]
    if "PassMethod" in kwds:
        argstrs.append(convert(kwds.pop("PassMethod")))
    argstrs += [", ".join((repr(k), convert(v))) for k, v in kwds.items()]
    return "{:>15}({});...".format(m_name(elem.__class__), ", ".join(argstrs))


def save_m(ring: Lattice, filename: str | None = None) -> None:
    """Save a :py:class:`.Lattice` as a Matlab m-file

    Parameters:
        ring:           Lattice description
        filename:       Name of the '.m' file. Default: outputs on
          :py:obj:`sys.stdout`

    See Also:
        :py:func:`.save_lattice` for a generic lattice-saving function.
    """

    def save(file):
        with np.printoptions(linewidth=1000, floatmode="unique"):
            print("ring = {...", file=file)
            for elem in matlab_ring(ring):
                print(_element_to_m(elem), file=file)
            print("};", file=file)

    if filename is None:
        save(sys.stdout)
    else:
        with open(filename, "w") as mfile:
            [funcname, _] = splitext(basename(filename))
            print(f"function ring = {funcname}()", file=mfile)
            save(mfile)
            print("    function v=False()\n        v=false;\n    end", file=mfile)
            print("    function v=True()\n        v=true;\n    end", file=mfile)
            print("end", file=mfile)


# Simulates the deprecated "mat_file" and "mat_key" attributes
def _mat_file(ring):
    """.mat input file. Deprecated, use *in_file* instead."""
    try:
        in_file = ring.in_file
    except AttributeError:
        raise AttributeError("'Lattice' object has no attribute 'mat_file'") from None
    if isinstance(in_file, str):
        _, ext = os.path.splitext(in_file)
        if ext != ".mat":
            raise AttributeError("'Lattice' object has no attribute 'mat_file'")
    else:
        raise AttributeError("'Lattice' object has no attribute 'mat_file'")
    return in_file


def _mat_key(ring):
    """selected Matlab variable. Deprecated, use *use* instead."""
    try:
        mat_key = ring.use
    except AttributeError:
        raise AttributeError("'Lattice' object has no attribute 'mat_key'") from None
    return mat_key


# noinspection PyUnusedLocal
def _ignore(ring, value):
    pass


register_format(
    ".mat",
    load_mat,
    save_mat,
    descr="Matlab binary mat-file. See :py:func:`.load_mat`.",
)

register_format(
    ".m",
    load_m,
    save_m,
    descr="Matlab text m-file. See :py:func:`.load_m`.",
)

Lattice.mat_file = property(_mat_file, _ignore, None)
Lattice.mat_key = property(_mat_key, _ignore, None)
