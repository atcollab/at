"""Lattice object

The methods implemented in this module are internal to the 'lattice' package.
This is necessary to ensure that the 'lattice' package is independent of other
AT packages.

Other Lattice methods are implemented in other AT packages and are available
as soon as the package is imported. The 'tracking' and 'physics' packages are
automatically imported.

As an example, see the at.physics.orbit module
"""

from __future__ import annotations

__all__ = [
    "Lattice",
    "Filter",
    "type_filter",
    "params_filter",
    "lattice_filter",
    "elem_generator",
    "no_filter",
]

import sys
import copy
import math

if sys.version_info.minor < 9:
    from typing import Callable, Iterable, Generator

    SupportsIndex = int
else:
    from typing import SupportsIndex
    from collections.abc import Callable, Iterable, Generator
import warnings

import numpy as np

from . import elements as elt
from .elements import Element
from .particle_object import Particle
from .utils import AtError, AtWarning, Refpts
from .utils import get_s_pos, get_elements, get_value_refpts, set_value_refpts

# noinspection PyProtectedMember
from .utils import get_uint32_index, get_bool_index, _refcount, Uint32Refpts
from .utils import refpts_iterator, checktype, get_geometry
from .transformation import set_rotation, set_tilt, set_shift
from ..constants import clight, e_mass

_TWO_PI_ERROR = 1.0e-4
Filter = Callable[..., Iterable[Element]]

_DEFAULT_PASS = {
    False: (
        ("cavity_pass", elt.RFCavity, "auto"),
        ("dipole_pass", elt.Dipole, "auto"),
        ("quadrupole_pass", elt.Quadrupole, "auto"),
        ("wiggler_pass", elt.Wiggler, "auto"),
        ("sextupole_pass", elt.Sextupole, "auto"),
        ("octupole_pass", elt.Octupole, "auto"),
        ("multipole_pass", elt.Multipole, "auto"),
        ("collective_pass", elt.Collective, "auto"),
        ("diffusion_pass", elt.QuantumDiffusion, "auto"),
        ("energyloss_pass", elt.EnergyLoss, "auto"),
        ("simplequantdiff_pass", elt.SimpleQuantDiff, "auto"),
        ("simpleradiation_pass", elt.SimpleRadiation, "auto"),
    ),
    True: (
        ("cavity_pass", elt.RFCavity, "auto"),
        ("dipole_pass", elt.Dipole, "auto"),
        ("quadrupole_pass", elt.Quadrupole, "auto"),
        ("wiggler_pass", elt.Wiggler, "auto"),
        ("sextupole_pass", elt.Sextupole, None),
        ("octupole_pass", elt.Octupole, None),
        ("multipole_pass", elt.Multipole, None),
        ("collective_pass", elt.Collective, "auto"),
        ("diffusion_pass", elt.QuantumDiffusion, "auto"),
        ("energyloss_pass", elt.EnergyLoss, "auto"),
        ("simplequantdiff_pass", elt.SimpleQuantDiff, "auto"),
        ("simpleradiation_pass", elt.SimpleRadiation, "auto"),
    ),
}

# Don't warn on floating-point errors
np.seterr(divide="ignore", invalid="ignore")
warnings.filterwarnings("always", category=AtWarning, module=__name__)


# noinspection PyAttributeOutsideInit
class Lattice(list):
    """Lattice object.

    An AT lattice is a sequence of AT elements. In addition to :py:class:`list`
    methods, :py:class:`Lattice` supports
    `extended indexing
    <https://numpy.org/doc/stable/user/basics.indexing.html#advanced-indexing>`_
    as a numpy :py:obj:`~numpy.ndarray`.
    """

    # Attributes displayed:
    _disp_attributes = (
        "name",
        "energy",
        "particle",
        "periodicity",
        "harmonic_number",
        "beam_current",
        "nbunch",
    )
    # excluded attributes
    _excluded_attributes = ("nbunch",)
    # Attributes propagated in copies:
    _std_attributes = (
        "name",
        "_energy",
        "_particle",
        "periodicity",
        "_cell_harmnumber",
        "_radiation",
        "_beam_current",
        "_fillpattern",
    )

    # noinspection PyUnusedLocal
    def __init__(self, *args, iterator: Filter = None, scan: bool = False, **kwargs):
        # noinspection PyUnresolvedReferences
        """
        Lattice(elements: Iterable[Element], **params)
        Lattice(filter, [filter, ...,] iterator=iter, **params)

        Parameters:
            elements: iterable of Element objects.
            iter: function called as :pycode:`iter(params, *args)`. It must
              return an iterable of :py:class:`.Element` objects for building
              the lattice. It must also fill the *params* dictionary providing
              the Lattice attributes.
            params: dictionary of lattice parameters. A custom iterator may
              add, remove or modify parameters. Finally, the remaining
              parameters will be set as Lattice attributes.

        Keyword Arguments:
            name (str): Name of the lattice, Default: "".
            energy (float):  Energy of the beam.
            periodicity (int): Number of periods. Default: 1. If <= 0, it will be
              deduced from the sum of bending angles.
            particle (Particle | str): circulating particle. May be "relativistic",
              "electron", "positron", "proton", "antiproton", "posmuon", "negmuon"
              or a Particle object. Default: "relativistic".
            iterator:  custom iterator (see below). Default :py:obj:`None`.
            *: all other keywords will be set as attributes of the Lattice object.
            beam_current:   Total current in the beam, used for collective effects [A].

        An iterator ``it`` is called as :pycode:`it(params, *args)` where ``args``
        and ``params`` are the arguments of the ``Lattice`` constructor. It
        must yield the AT ``Elements`` for building the lattice. It must
        also fill its ``params`` dictionary argument, which will be used to
        set the ``Lattice`` attributes.
        An iterator can be:

        - a "generator" which yields elements from scratch.
          Examples: a list, or a file iterator,
        - a "filter" which runs through an input iterator, processes each
          element, possibly adds parameters to the params dictionary
          and yields the processed elements.

        .. Note::

           To reduce the inter-package dependencies, some methods of the
           lattice object are defined in other AT packages, in the module where
           the underlying function is implemented.

        .. Note::

           It is possible to define a filling pattern for the beam using the
           function :pycode:`ring.set_fillingpattern()`. The default configuration
           (no arguments) is for single bunch and is the one loaded at lattice
           initialization. See function help for details.
           Changing ``Lattice.harmonic_number`` will reset the filling pattern
           to its default configuration.
           The beam current can be changed with :pycode:`Lattice.beam_current = current`
           The filling pattern and beam current are used by collective effects
           passmethods.


        Examples:
            Chaining iterators (taken from ``load_mat``):

            >>> ring = Lattice(ringparam_filter, matfile_generator, filename,
                               iterator=params_filter, **params)

            ``matfile_generator(params, filename)``
                opens filename and generates AT elements for each cell of the
                Matlab cell array representing the lattice,

            ``ringparam_filter(params, matfile_generator, *args)``
                runs through ``matfile_generator(params, *args)``, looks for
                RingParam elements, fills params with their information and
                discards them,

            ``params_filter(params, ringparam_filter, *args)``
                runs through ``ringparam_filter(params, *args)``, looks for
                energy and periodicity if not yet defined.
        """
        if iterator is None:
            (arg1,) = args or [[]]  # accept 0 or 1 argument
            if isinstance(arg1, Lattice):
                elems = lattice_filter(kwargs, arg1)
            else:
                elems = params_filter(kwargs, type_filter, arg1)
        else:
            elems = iterator(kwargs, *args)

        super().__init__(elems)

        # removing excluded attributes
        for attr in self._excluded_attributes:
            kwargs.pop(attr, None)
        # set default values
        kwargs.setdefault("name", "")
        periodicity = kwargs.setdefault("periodicity", 1)
        kwargs.setdefault("particle", kwargs.pop("_particle", Particle()))
        kwargs.setdefault("beam_current", kwargs.pop("_beam_current", 0.0))
        # dummy initialization in case the harmonic number is not there
        kwargs.setdefault("_fillpattern", np.ones(1))
        # Remove temporary keywords
        frequency: float | None = kwargs.pop("_frequency", None)
        cell_length: float | None = kwargs.pop("_length", None)
        cell_h = kwargs.pop("cell_harmnumber", kwargs.pop("_cell_harmnumber", math.nan))
        ring_h = kwargs.pop("harmonic_number", cell_h * periodicity)

        energy = kwargs.setdefault("energy", kwargs.pop("_energy", None))
        if energy is None:
            raise AtError("Lattice energy is not defined")

        # set attributes
        self.update(kwargs)

        # Setting the harmonic number is delayed to have self.beta available
        if not (frequency is None or frequency == 0.0):
            rev = self.beta * clight / cell_length
            self._cell_harmnumber = int(round(frequency / rev))
            try:
                fp = kwargs.pop("_fillpattern", np.ones(1))
                self.set_fillpattern(bunches=fp)
            except AssertionError:
                self.set_fillpattern()
        elif not ((ring_h is None) or math.isnan(ring_h)):
            self._cell_harmnumber = ring_h / periodicity

    def __getitem__(self, key):
        try:  # Integer
            return super().__getitem__(key.__index__())
        except (AttributeError, TypeError):
            if isinstance(key, slice):  # Slice
                rg = range(*key.indices(len(self)))
            else:  # Array of integers or boolean
                rg = get_uint32_index(self, key, endpoint=False)
            return Lattice(
                elem_generator,
                (super(Lattice, self).__getitem__(i) for i in rg),
                iterator=self.attrs_filter,
            )

    def __setitem__(self, key, values):
        try:  # Integer or slice
            super().__setitem__(key, values)
        except TypeError:  # Array of integers or boolean
            rg = get_uint32_index(self, key, endpoint=False)
            for i, v in zip(*np.broadcast_arrays(rg, values)):
                super().__setitem__(i, v)

    def __delitem__(self, key):
        try:  # Integer or slice
            super().__delitem__(key)
        except TypeError:  # Array of integers or boolean
            rg = get_uint32_index(self, key, endpoint=False)
            for i in reversed(rg):
                super().__delitem__(i)

    def __repr__(self):
        ks = ", ".join(f"{k}={v!r}" for k, v in self.attrs.items())
        return f"Lattice({super().__repr__()}, {ks})"

    def __str__(self):
        ks = ", ".join(f"{k}={v!r}" for k, v in self.attrs.items())
        return f"Lattice(<{len(self)} elements>, {ks})"

    def __add__(self, elems):
        """Add elems, an iterable of AT elements, to the lattice"""
        return self.concatenate(elems, copy=True)

    def __iadd__(self, elems):
        self.concatenate(elems, copy=False)
        return self

    def __mul__(self, n):
        return self.repeat(n)

    def _addition_filter(self, elems: Iterable[Element], copy_elements=False):
        cavities = []
        length = 0.0
        params = {}

        for elem in type_filter(params, elems):
            if isinstance(elem, elt.RFCavity):
                cavities.append(elem)
                elem.Energy = self._energy
            elif elem.PassMethod.endswith("RadPass"):
                elem.Energy = self._energy
            elif hasattr(elem, "Energy"):
                del elem.Energy
            elif hasattr(elem, "_turnhistory"):
                elem.clear_history(self)
            length += getattr(elem, "Length", 0.0)
            if copy_elements:
                yield elem.deepcopy()
            else:
                yield elem

        if cavities and not hasattr(self, "_cell_harmnumber"):
            cavities.sort(key=lambda el: el.Frequency)
            try:
                self._cell_harmnumber = cavities[0].HarmNumber
            except AttributeError:
                length += self.get_s_pos(len(self))[0]
                rev = self.beta * clight / length
                frequency = cavities[0].Frequency
                self._cell_harmnumber = int(round(frequency / rev))
        self._radiation |= params.pop("_radiation")

    def insert(self, idx: SupportsIndex, elem: Element, copy_elements=False):
        r"""This method allow to insert an AT element in the lattice.

        Parameters:
            idx (SupportsIndex): index at which the lement is inserted
            elem (Element): AT element to be inserted in the lattice
            copy_elements(bool): Default :py:obj:`True`.
                                 If :py:obj:`True` a deep copy of elem
                                 is used.
        """
        # scan the new element to update it
        elist = list(  # noqa: F841
            self._addition_filter([elem], copy_elements=copy_elements)
        )
        super().insert(idx, elem)

    def extend(self, elems: Iterable[Element], copy_elements=False):
        # noinspection PyUnresolvedReferences
        r"""This method adds all the elements of `elems` to the end of the
        lattice. The behavior is the same as for a :py:obj:`list`

        Equivalent syntaxes:

        >>> ring.extend(elems)
        >>> ring += elems

        Parameters:
            elems (Iterable[Element]): Sequence of AT elements to be
                                      appended to the lattice
            copy_elements(bool): Default :py:obj:`True`.
                                 If :py:obj:`True` deep copies of each
                                 element of elems are used
        """
        if hasattr(self, "_energy"):
            # When unpickling a Lattice, extend is called before the lattice
            # is initialized. So skip this.
            elems = self._addition_filter(elems, copy_elements=copy_elements)
        super().extend(elems)

    def append(self, elem: Element, copy_elements=False):
        # noinspection PyUnresolvedReferences
        r"""This method overwrites the inherited method :py:meth:`list.append()`,
        its behavior is changed, it accepts only AT lattice elements
        :py:obj:`Element` as input argument.

        Equivalent syntaxes:

        >>> ring.append(elem)
        >>> ring += [elem]

        Parameters:
            elem (Element): AT element to be appended to the lattice
            copy_elements(bool): Default :py:obj:`True`.
                                 If :py:obj:`True` a deep copy of elem
                                 is used
        """
        self.extend([elem], copy_elements=copy_elements)

    def repeat(self, n: int, copy_elements: bool = True):
        # noinspection SpellCheckingInspection,PyUnresolvedReferences,PyRedeclaration
        r"""This method allows to repeat the lattice `n` times.
        If `n` does not divide `ring.periodicity`, the new ring
        periodicity is set to 1, otherwise  it is set to
        `ring.periodicity /= n`.

        Equivalent syntaxes:

        >>> newring = ring.repeat(n)
        >>> newring = ring * n

        Parameters:
            n : number of repetitions
            copy_elements:  If :py:obj:`True`, deep copies of the lattice are used for
              the repetition. Otherwise, the original elements are repeated in the
              developed lattice.

        Returns:
            newring (Lattice): the new repeated lattice
        """

        def copy_fun(elem, copy):
            if copy:
                return elem.deepcopy()
            else:
                return elem

        periodicity = self.periodicity
        if n != 0 and periodicity > 1:
            nbp = periodicity / n
            periodicity = int(round(nbp))
            if abs(periodicity - nbp) > _TWO_PI_ERROR:
                warnings.warn(
                    AtWarning(
                        f"Non-integer number of cells: {self.periodicity}/{n}."
                        " Periodicity set to 1"
                    ),
                    stacklevel=1,
                )
                periodicity = 1
        hdict = {"periodicity": periodicity}
        try:
            hdict.update(harmonic_number=self.cell_harmnumber * n * periodicity)
        except AttributeError:
            pass
        elems = (copy_fun(el, copy_elements) for _ in range(n) for el in self)
        return Lattice(elem_generator, elems, iterator=self.attrs_filter, **hdict)

    def concatenate(
        self, *lattices: Iterable[Element], copy_elements=False, copy=False
    ):
        # noinspection PyUnresolvedReferences,SpellCheckingInspection,PyRedeclaration
        """Concatenate several `Iterable[Element]` with the lattice

        Equivalent syntaxes:

        >>> newring = ring.concatenate(r1, r2, r3, copy=True)
        >>> newring = ring + r1 + r2 + r3

        >>> ring.concatenate(r1, copy=False)
        >>> ring += r1

        Parameters:
            lattices: :py:obj:`Iterables[Element]` to be concatenanted
                      to the Lattice, several lattices are allowed
                      (see example)
            copy_elements(bool): Default :py:obj:`False`. If :py:obj:`True`
                        deepcopies of the elements of lattices are used
            copy(bool): Default :py:obj:`False`. If :py:obj:`True`
                           the lattice is modified in place.
                           Oterwise a new Lattice object is returned

        Returns:
            lattice(Lattice): concatenated Lattice, if `copy==True` the
                              new lattice object is returned
                              otherwise None
        """
        if copy:
            lattice = Lattice(self)
        else:
            lattice = self
        for lat in lattices:
            lattice.extend(lat, copy_elements=copy_elements)
        return lattice if copy else None

    def reverse(self, copy=False):
        # noinspection PyUnresolvedReferences
        r"""Reverse the order of the lattice and swapt the faces
        of elements. Alignment errors are not swapped

        Parameters:
            copy(bool): Default :py:obj:`False`. If :py:obj:`True`
                           the lattice is modified in place.
                           Oterwise a new Lattice object is returned

        Returns:
            lattice(Lattice): reversed Lattice, if `copy==True` the
                              new lattice object is returned
                              otherwise None
        Example:
            >>> newring = ring.reverse(copy=True)
        """
        elems = (el.swap_faces(copy=True) for el in reversed(self))
        if copy:
            return Lattice(elem_generator, elems, iterator=self.attrs_filter)
        else:
            reversed_list = list(elems)
            self[:] = reversed_list

    def develop(self, copy_elements: bool = True) -> Lattice:
        """Develop a periodical lattice by repeating its elements
        *self.periodicity* times

        Parameters:
            copy_elements:  If :py:obj:`True`, deep copies of the elements are used for
              the repetition. Otherwise, the original elements are repeated in the
              developed lattice.

        Returns:
            newlattice: The developed lattice
        """
        return self.repeat(self.periodicity, copy_elements=copy_elements)

    @property
    def attrs(self) -> dict:
        """Dictionary of lattice attributes"""

        def extattr(d):
            for k in self._disp_attributes:
                d.pop(k, None)
                yield k, getattr(self, k, None)

        vrs = vars(self).copy()
        # Standard attributes
        res = {k: v for k, v in extattr(vrs) if v is not None}
        # Custom attributes
        res.update((k, v) for k, v in vrs.items() if not k.startswith("_"))
        return res

    def rotate(self, n: int) -> Lattice:
        """Return a new lattice rotated left by n elements"""
        if len(self) == 0:
            return self.copy()
        n = n % len(self)  # works for n<0
        return self[n:] + self[:n]

    def update(self, *args, **kwargs) -> None:
        """
        update(**kwargs)
        update(mapping, **kwargs)
        update(iterable, **kwargs)

        Update the lattice attributes with the given values
        """
        attrs = dict(*args, **kwargs)
        for key, value in attrs.items():
            setattr(self, key, value)

    def copy(self) -> Lattice:
        """Returns a shallow copy of the lattice"""
        return copy.copy(self)

    def deepcopy(self) -> Lattice:
        """Returns a deep copy of the lattice"""
        return copy.deepcopy(self)

    def slice_elements(self, refpts: Refpts, slices: int = 1) -> Lattice:
        """Create a new lattice by slicing the elements at refpts

        Parameters:
            refpts:     Element selector
            slices:     Number of slices in the specified range. Ignored if
              size is specified. Default: no slicing

        Returns:
            newring:    New Lattice object
        """

        def slice_generator(_):
            check = get_bool_index(self, refpts)
            for el, ok in zip(self, check):
                if ok and (slices > 1):
                    frac = np.ones(slices) / slices
                    yield from el.divide(frac)
                else:
                    yield el

        return Lattice(slice_generator, iterator=self.attrs_filter)

    def slice(self, size: float | None = None, slices: int | None = 1) -> Lattice:
        """Create a new lattice by slicing the range of interest into small
        elements

        Keyword arguments:
            size=None:      Length of a slice. Default: computed from the
              range and number of points: ``size = (s_max-s_min)/slices``.
            slices=1:       Number of slices in the specified range. Ignored if
              size is specified. Default: no slicing

        Returns:
            newring:    New Lattice object
        """

        def slice_generator(_, ibeg, iend):
            yield from self[:ibeg]
            for elem in self[ibeg:iend]:
                nslices = int(math.ceil(elem.Length / size))
                if nslices > 1:
                    frac = np.ones(nslices) / nslices
                    yield from elem.divide(frac)
                else:
                    yield elem
            yield from self[iend:]

        s_range = self.s_range
        if size is None:
            smin, smax = s_range
            size = (smax - smin) / slices
        i_range = self.i_range
        return Lattice(
            slice_generator,
            i_range[0],
            i_range[-1],
            iterator=self.attrs_filter,
            s_range=s_range,
        )

    def attrs_filter(
        self, params: dict, elem_filter: Filter, *args
    ) -> Iterable[Element]:
        """Filter function which duplicates the lattice attributes

        Parameters:
            params:         Dictionary of Lattice attributes
            elem_filter:    Element filter
            *args:          Arguments provided to *elem_filter*
        """
        for key in self._std_attributes:
            try:
                params.setdefault(key, getattr(self, key))
            except AttributeError:
                pass
        return elem_filter(params, *args)

    @property
    def s_range(self) -> None | tuple[float, float]:
        """Range of interest: (s_min, s_max). :py:obj:`None` means the full cell."""
        try:
            return self._s_range
        except AttributeError:
            self.s_range = None
            return self._s_range

    # noinspection PyAttributeOutsideInit
    @s_range.setter
    def s_range(self, value: tuple[float, float] | None):
        spos = self.get_s_pos(range(len(self) + 1))
        if value is None:
            value = (0.0, spos[-1])
        ok = np.flatnonzero(np.logical_and(spos > value[0], spos < value[1]))
        if len(ok) > 0:
            i1 = max(ok[0] - 1, 0)
            i2 = min(ok[-1] + 1, len(self))
            self._i_range = range(i1, i2 + 1)
        else:
            self._i_range = range(0, 0)
        self._s_range = value

    @property
    def i_range(self) -> Uint32Refpts:
        """Range of elements inside the range of interest"""
        try:
            i_range = self._i_range
        except AttributeError:
            self.s_range = None
            i_range = self._i_range
        return get_uint32_index(self, i_range)

    @property
    def energy(self) -> float:
        """Lattice energy"""
        return self._energy

    @energy.setter
    def energy(self, energy: float):
        # Set the Energy attribute of radiating elements
        for elem in self:
            if isinstance(
                elem, (elt.RFCavity, elt.Wiggler)
            ) or elem.PassMethod.endswith("RadPass"):
                elem.Energy = energy
        # Set the energy attribute of the Lattice
        # Use a numpy scalar to allow division by zero
        # self._energy = np.array(energy, dtype=float)
        self._energy = energy

    @property
    def cell_length(self) -> float:
        """Cell length (1 cell) [m]

        See Also:

            :py:meth:`circumference`.
        """
        return self.get_s_pos(len(self))[0]

    @property
    def cell_revolution_frequency(self) -> float:
        """Revolution frequency for 1 cell [Hz]

        See Also:

            :py:meth:`revolution_frequency`.
        """
        beta = self.beta
        return beta * clight / self.cell_length

    @property
    def cell_harmnumber(self) -> float:
        """Harmonic number per cell

        See Also:

            :py:meth:`harmonic_number`.
        """
        try:
            return self._cell_harmnumber
        except AttributeError as exc:
            raise AttributeError(
                "harmonic_number undefined: No cavity found in the lattice"
            ) from exc

    @property
    def circumference(self) -> float:
        """Ring circumference (full ring) [m]

        See Also:

            :py:meth:`cell_length`.
        """
        return self.periodicity * self.get_s_pos(len(self))[0]

    @property
    def revolution_frequency(self) -> float:
        """Revolution frequency (full ring) [Hz]

        See Also:

            :py:meth:`cell_revolution_frequency`.
        """
        beta = self.beta
        return beta * clight / self.circumference

    @property
    def particle(self) -> Particle:
        """Circulating particle"""
        return self._particle

    @particle.setter
    def particle(self, particle: str | Particle):
        if isinstance(particle, Particle):
            self._particle = particle
        else:
            self._particle = Particle(particle)

    def set_wake_turnhistory(self):
        """Function to reset the shape of the turn history
        in collective elements based on the number of slices,
        turns and bunches
        """
        for e in self:
            if e.is_collective:
                e.clear_history(ring=self)

    def set_fillpattern(self, bunches: int | np.ndarray = 1):
        """Function to generate the filling pattern lof the ring.
        The filling pattern is computed as:

        ``bunches/numpy.sum(bunches)``

        This function also generates the bunch spatial distribution
        accessible with ``Lattice.bunch_spos``

        Keyword Arguments:
           bunches:  integer or array of positive double or bool to define
                     the bunch distribution.
                     For scalar input, equidistant bunches are assumed.
                     ``ring.harmonic_number`` has to be a multiple of
                     ``bunches``.
                     For array input the condition
                     ``len(bunches)==ring.harmonic_number`` is required.
                     (default=1, single bunch configuration).
        """
        if isinstance(bunches, int):
            if self.harmonic_number % bunches == 0:
                fp = np.zeros(self.harmonic_number)
                fp[:: int(self.harmonic_number / bunches)] = 1
            else:
                raise AtError(
                    "Harmonic number has to be a multiple of the scalar input bunches"
                )
        elif np.isscalar(bunches):
            raise AtError("Scalar input for bunches must be an integer")
        else:
            bunches = bunches.astype(dtype=float, casting="safe", copy=False)
            assert (
                len(bunches) == self.harmonic_number
            ), f"bunches array input has to be of shape ({self.harmonic_number},)"
            assert np.all(
                bunches >= 0.0
            ), "bunches array can contain only positive numbers"
            fp = bunches
        self._fillpattern = fp / np.sum(fp)
        self.set_wake_turnhistory()

    @property
    def fillpattern(self) -> np.ndarray:
        """Filling pattern describing the bunch relative
        amplitudes such that ``sum(fillpattern)=1``
        """
        return self._fillpattern

    @fillpattern.setter
    def fillpattern(self, value):
        """Filling pattern describing the bunch relative
        amplitudes such that ``sum(fillpattern)=1``.
        Calls the function ``Lattice.set_fillpattern``.
        """
        self.set_fillpattern(value)

    @property
    def bunch_list(self) -> np.ndarray:
        """Indices of filled bunches"""
        return np.flatnonzero(self._fillpattern)

    def get_beam_current(self):
        return self._beam_current

    def set_beam_current(self, value, clear_history=True):
        self._beam_current = value
        if clear_history:
            self.set_wake_turnhistory()

    beam_current = property(get_beam_current, set_beam_current)

    @property
    def bunch_currents(self) -> np.ndarray:
        """Bunch currents [A]"""
        return self.beam_current * self._fillpattern[self._fillpattern > 0]

    @property
    def bunch_spos(self) -> np.ndarray:
        """Bunch position around the ring [m]"""
        try:
            circ = self.beta * clight * self.harmonic_number / self.rf_frequency
        except AtError:
            circ = self.circumference
        bs = circ / len(self._fillpattern)
        allpos = bs * np.arange(len(self._fillpattern))
        return allpos[self._fillpattern > 0]

    @property
    def nbunch(self) -> int:
        """Number of bunches"""
        return np.count_nonzero(self._fillpattern)

    @property
    def harmonic_number(self) -> int:
        """Ring harmonic number (full ring)

        See Also:

            :py:meth:`cell_harmnumber`.
        """
        try:
            return int(self.periodicity * self._cell_harmnumber)
        except AttributeError as exc:
            raise AttributeError(
                "harmonic_number undefined: " "No cavity found in the lattice"
            ) from exc

    @harmonic_number.setter
    def harmonic_number(self, value):
        cell_h = float(value) / self.periodicity
        # check on ring
        if value - round(value) != 0:
            raise AtError(f"harmonic number ({value}) must be integer")
        # check on cell
        # if cell_h-round(cell_h) != 0:
        #     raise AtError('harmonic number ({}) must be a multiple of {}'
        #                   .format(value, int(self.periodicity)))
        self._cell_harmnumber = cell_h
        if len(self._fillpattern) != value:
            warnings.warn(
                AtWarning(
                    "Harmonic number changed, resetting fillpattern to "
                    "default (single bunch)."
                ),
                stacklevel=1,
            )
            self.set_fillpattern()

    @property
    def gamma(self) -> float:
        r"""Relativistic :math:`\gamma` of the particles"""
        rest_energy = self.particle.rest_energy
        if rest_energy == 0.0:
            rest_energy = e_mass
        return float(self.energy / rest_energy)

    @property
    def beta(self) -> float:
        r"""Relativistic :math:`\beta` of the particles"""
        rest_energy = self.particle.rest_energy
        if rest_energy == 0.0:
            return 1.0
        else:
            gamma = float(self.energy / rest_energy)
            return math.sqrt(1.0 - 1.0 / gamma / gamma)

    # noinspection PyPep8Naming
    @property
    def BRho(self) -> float:
        """Magnetic rigidity [T.m]"""
        rest_energy = self.particle.rest_energy
        if rest_energy == 0.0:
            rest_energy = e_mass
        return math.sqrt(self.energy**2 - rest_energy**2) / clight

    @property
    def is_6d(self) -> bool:
        """True if at least one element modifies the beam momentum

        See Also:

            :py:meth:`enable_6d`, :py:meth:`disable_6d`.
        """
        try:
            return self._radiation
        except AttributeError:
            radiate = False
            for elem in self:
                if elem.longt_motion:
                    radiate = True
                    break
            # noinspection PyAttributeOutsideInit
            self._radiation = radiate
            return radiate

    @property
    def is_collective(self) -> bool:
        """:py:obj:`True` if any element involves collective effects"""
        for elem in self:
            if elem.is_collective:
                return True
        return False

    @property
    def has_cavity(self) -> bool:
        """:py:obj:`True` if the lattice contains an active
        :py:class:`RFCavity`
        """
        for elem in self:
            if elem.PassMethod.endswith("CavityPass"):
                return True
        return False

    # noinspection PyShadowingNames
    def modify_elements(
        self, elem_modify: Callable, copy: bool | None = True, **kwargs
    ):
        """Modify selected elements, in-place or in a lattice copy

        Parameters:
            elem_modify : element selection function.
              If ``elem_modify(elem)`` returns :py:obj:`None`, the element is
              unchanged. Otherwise, ``elem_modify(elem)`` must return a
              dictionary of attribute name and values, to be set to elem.
            copy:         If True, return a shallow copy of the lattice.
              Only the modified elements are copied.
              If False, the modification is done in-place

        Returns:
            newring:    New lattice if copy is :py:obj:`True`, :py:obj:`None`
              if copy is :py:obj:`False`

        Keyword Arguments:
        """

        def lattice_modify(**kws):
            """Modifies the Lattice with elements modified by elem_modify"""
            radiate = False
            for elem in self:
                attrs = elem_modify(elem)
                if attrs is not None:
                    elem.update(attrs)
                if elem.longt_motion:
                    radiate = True
            self._radiation = radiate
            self.update(kws)

        def lattice_copy(params):
            """Custom iterator for the creation of a new lattice"""
            radiate = False
            for elem in self:
                attrs = elem_modify(elem)
                if attrs is not None:
                    elem = elem.copy()
                    elem.update(attrs)
                if elem.longt_motion:
                    radiate = True
                yield elem
            params["_radiation"] = radiate

        if copy:
            return Lattice(lattice_copy, iterator=self.attrs_filter, **kwargs)
        else:
            lattice_modify(**kwargs)

    def _set_6d(self, enable: bool, *args, **kwargs):
        """Set the lattice radiation state"""

        def lattice_modify():
            """Modifies the Lattice in place"""
            radiate = False
            for elem in self:
                new_pass = getpass(elem)
                if new_pass:
                    elem.set_longt_motion(enable, new_pass=new_pass, **vargs)
                if elem.longt_motion:
                    radiate = True
            self._radiation = radiate
            self.update(kwargs)

        def lattice_copy(params):
            """Custom iterator for the creation of a new lattice"""
            radiate = False
            for elem in self:
                new_pass = getpass(elem)
                if new_pass:
                    elem = elem.set_longt_motion(
                        enable, new_pass=new_pass, copy=True, **vargs
                    )
                if elem.longt_motion:
                    radiate = True
                yield elem
            params["_radiation"] = radiate

        cp = kwargs.pop("copy", False)
        if len(args) > 0:

            def getpass(elem):
                return "auto" if isinstance(elem, args) else None

            if not all(issubclass(cl, elt.LongtMotion) for cl in args):
                raise TypeError("All arguments must be subclasses of 'LongtMotion'")
            if len(kwargs) > 0:
                raise AtError("No keyword is allowed in this mode")
        else:

            def getpass(elem):
                for eltype, psm in pass_table:
                    if isinstance(elem, eltype):
                        return psm
                return None

            def passm(key, eltype, def_pass):
                if allset:
                    def_pass = all_pass
                return eltype, kwargs.pop(key, def_pass)

            # Look for global defaults
            try:
                all_pass = kwargs.pop("all_pass")
            except KeyError:
                allset = False
            else:
                allset = True
            # Build table of PassMethods
            pass_table = [passm(*defs) for defs in _DEFAULT_PASS[enable]]

        vargs = {"energy": self.energy} if enable else {}
        if cp:
            return Lattice(lattice_copy, iterator=self.attrs_filter, **kwargs)
        else:
            lattice_modify()

    # noinspection PyShadowingNames,PyIncorrectDocstring
    def enable_6d(self, *args, **kwargs) -> Lattice | None:
        # noinspection PyUnresolvedReferences
        r"""
        enable_6d(elem_class[, elem_class]..., copy=False)
        enable_6d(cavity_pass='auto'[, dipole_pass='auto']..., copy=False)

        Turn longitudinal motion on. By default, turn on
        radiation in dipoles and quadrupoles, turn on RF cavities, activates
        collective effects and other elements acting on momentum.

        Modify the lattice in-place or creates a new lattice, depending on the
        ``copy`` keyword argument.

        **Syntax using positional arguments:**

        Parameters:
            elem_class:                 :py:class:`.LongtMotion` subclass.
              Longitudinal motion is turned on for elements which are
              instances of any given ``elem_class``.

              In adition to single element classes, a few grouping classes are
              available:

              * :py:class:`.LongtMotion`: all elements possibly acting on
                momentum,
              * :py:class:`.Radiative`: default radiative elements:
                :py:class:`.Dipole`, :py:class:`.Quadrupole`,
                :py:class:`.Wiggler`,
              * :py:class:`.Collective`: all elements dealing with collective
                effects.

              The default PassMethod conversion is used, as with the ``'auto'``
              keyword value..

              **No keyword except** ``copy`` **is allowed in this syntax.**

        **Syntax using keyword arguments:**

        Keyword arguments:
            all_pass:                   PassMethod overloading the default
              values for all elements (:py:obj:`None` or 'auto')
            cavity_pass='auto':         PassMethod set on cavities
            dipole_pass='auto':         PassMethod set on dipoles
            quadrupole_pass='auto':     PassMethod set on quadrupoles
            wiggler_pass='auto':        PassMethod set on wigglers
            sextupole_pass=None:        PassMethod set on sextupoles
            octupole_pass=None:         PassMethod set on octupoles
            multipole_pass=None  :      PassMethod set on higher order
              multipoles
            collective_pass='auto':     PassMethod set on collective effect
              elements (:py:class:`.WakeElement`,...)
            diffusion_pass='auto':      PassMethod set on
              :py:class:`.QuantumDiffusion`
            copy=False: If :py:obj:`False`, the modification is done in-place,
              If ``True``, return a shallow copy of the lattice. Only the
              radiating elements are copied with PassMethod modified.

              .. Caution::

                 a shallow copy means that all non-modified elements are shared
                 with the original lattice. Any further modification will
                 affect both lattices.

        For PassMethod names, the convention is:

        * :py:obj:`None`:        No change
        * 'auto':          Use the default conversion (replace \*Pass by
          \*RadPass)
        * anything else:   set as the new PassMethod

        Examples:
            >>> ring.enable_6d()

            Modify `ring` in-place, turn cavities ON, turn synchrotron
            radiation ON in Dipoles and Quadrupoles, turn collective effects
            ON.

            >>> ring.enable_6d(at.RFCavity, at.Radiative)

            Modify `ring` in-place, turn cavities ON, turn synchrotron
            radiation ON in dipoles and quadupoles.

            >>> newring = ring.enable_6d(at.Collective, copy=True)

            Returns a new lattice with collective effects turned ON and nothing
            else changed

            >>> newring = ring.enable_6d(all_pass=None, collective_pass="auto", copy=True)

            Same as the previous example, using the keyword syntax.

        See Also:

            :py:meth:`disable_6d`, :py:attr:`is_6d`.
        """
        return self._set_6d(True, *args, **kwargs)

    # noinspection PyShadowingNames,PyIncorrectDocstring
    def disable_6d(self, *args, **kwargs) -> Lattice | None:
        # noinspection PyUnresolvedReferences
        r"""
        disable_6d(elem_class[, elem_class]... , copy=False)
        disable_6d(cavity_pass='auto'[, dipole_pass='auto']..., copy=False)

        Turn longitudinal motion off. By default, remove all longitudinal
        motion.

        Modify the lattice in-place or creates a new lattice, depending on the
        ``copy`` keyword argument.

        **Syntax using positional arguments:**

        Parameters:
            elem_class:                 :py:class:`.LongtMotion` subclass.
              Longitudinal motion is turned off for elements which are
              instances of any given ``elem_class``.

              In adition to single element classes, a few grouping classes are
              available:

              * :py:class:`.LongtMotion`: all elements possibly acting on
                momentum,
              * :py:class:`.Radiative`: default radiative elements:
                :py:class:`.Dipole`, :py:class:`.Quadrupole`,
                :py:class:`.Wiggler`,
              * :py:class:`.Collective`: all elements dealing with collective
                effects.

              The default PassMethod conversion is used, as with the ``'auto'``
              keyword value.

              **No keyword except** ``copy`` **is allowed in this syntax.**

        **Syntax using keyword arguments:**

        Keyword arguments:
            all_pass:                   PassMethod overloading the default
              values for all elements (:py:obj:`None` or 'auto')
            cavity_pass='auto':         PassMethod set on cavities
            dipole_pass='auto':         PassMethod set on dipoles
            quadrupole_pass=auto:       PassMethod set on quadrupoles
            wiggler_pass='auto':        PassMethod set on wigglers
            sextupole_pass='auto':      PassMethod set on sextupoles
            octupole_pass='auto':       PassMethod set on octupoles
            multipole_pass='auto':      PassMethod set on higher order
              multipoles
            collective_pass='auto':     PassMethod set on collective effect
              elements (:py:class:`.WakeElement`,...)
            diffusion_pass='auto':      PassMethod set on
              :py:class:`.QuantumDiffusion`
            copy=False: If :py:obj:`False`, the modification is done in-place,
              If ``True``, return a shallow copy of the lattice. Only the
              radiating elements are copied with PassMethod modified.

              .. Caution::

                 a shallow copy means that all non-modified elements are shared
                 with the original lattice. Any further modification will
                 affect both lattices.

        For PassMethod names, the convention is:

        * :py:obj:`None`:        no change
        * 'auto':          Use the default conversion (replace \*RadPass by
          \*Pass)
        * anything else:   set as the new PassMethod

        Examples:
            >>> ring.disable_6d()

            Modify `ring` in-place, turn OFF everything affecting the
            longitudinal momentum.

            >>> ring.disable_6d(at.RFCavity)

            Turn cavities OFF (nothing else modified).

            >>> ring.disable_6d(all_pass=None, cavity_pass="auto")

            Same as the previous example, but using the keyword syntax.

            >>> newring = ring.disable_6d(cavity_pass=None, copy=True)

            Return a new Lattice (shallow copy of `ring`) with everything
            turned OFF except RF cavities.

            >>> newring = ring.disable_6d(
            ...     all_pass=None, sextupole_pass="DriftPass", copy=True
            ... )

            Return a new Lattice (shallow copy of `ring`) with sextupoles
            turned into Drifts (turned off) and everything else unchangedd.

        See Also:

            :py:meth:`enable_6d`, :py:attr:`is_6d`.
        """
        return self._set_6d(False, *args, **kwargs)

    def sbreak(self, break_s, break_elems=None, **kwargs):
        """Insert elements at selected locations in the lattice

        Parameters:
            break_s:        location or array of locations of breakpoints
            break_elems:    elements to be inserted at breakpoints (array of
                            elements as long as break_s or single element
                            duplicated as necessary). Default: Marker('sbreak')
        Returns:
            newring:    A new lattice with new elements inserted at breakpoints
        """

        def sbreak_iterator(_, itmk):
            """Iterate over elements and breaks where necessary"""

            def next_mk():
                """Extract the next element to insert"""
                try:
                    return next(itmk)
                except StopIteration:
                    return sys.float_info.max, None

            s_end = 0.0
            # get the 1st insertion
            smk, mk = next_mk()
            # skip all insertions at negative break_s, if any
            while smk < s_end:
                smk, mk = next_mk()

            for elem in self:
                s_end += elem.Length
                # loop over all insertions within the element
                while smk < s_end:
                    frac = (s_end - smk) / elem.Length
                    if frac < 1.0:  # breakpoint is within the element
                        el0, elem = elem.divide([1.0 - frac, frac])
                        yield el0
                    yield mk
                    smk, mk = next_mk()
                yield elem

        # set default insertion
        if break_elems is None:
            break_elems = elt.Marker("sbreak")
        break_elems = np.reshape(break_elems, -1)
        # Check element lengths
        if not all(e.Length == 0 for e in break_elems):
            warnings.warn(
                AtWarning("Inserting elements with length!=0 may change the lattice"),
                stacklevel=2,
            )
        # broadcast break_s and break_elems to arrays of same size
        # and create an iterator over the elements to be inserted
        iter_mk = zip(*np.broadcast_arrays(break_s, break_elems))

        return Lattice(sbreak_iterator, iter_mk, iterator=self.attrs_filter, **kwargs)

    def reduce(self, **kwargs) -> Lattice:
        """Removes all elements with an ``IdentityPass`` PassMethod and merges
        compatible consecutive elements.

        The "reduced" lattice has the same optics as the original one, but
        fewer elements. Compatible elements must have the same ``PolynomA``,
        ``PolynomB`` and bending radius, so that the optics is preserved. But
        magnet misalignments are not preserved, so this method should be
        applied to lattices without errors.

        Keyword Args:
            keep (Refpts):      Keep the selected elements, even with
                ``IdentityPass`` PassMethod. Default: keep :py:class:`.Monitor`
                and :py:class:`.RFCavity` elements
        """

        def reduce_filter(_, itelem):
            try:
                elem = next(itelem).copy()
            except StopIteration:
                return
            for nxt in itelem:
                try:
                    elem.merge(nxt)
                except TypeError:
                    yield elem
                    elem = nxt.copy()
            yield elem

        kp = get_bool_index(self, lambda el: el.PassMethod != "IdentityPass")
        try:
            keep = get_bool_index(self, kwargs.pop("keep"))
        except KeyError:
            keep = get_bool_index(self, checktype((elt.Monitor, elt.RFCavity)))

        return Lattice(
            reduce_filter, self.select(kp | keep), iterator=self.attrs_filter, **kwargs
        )

    def replace(self, refpts: Refpts, **kwargs) -> Lattice:
        """Return a shallow copy of the lattice replacing the selected
        elements by a deep copy

        Parameters:
            refpts: element selector
        """
        check = get_bool_index(self, refpts)
        elems = (el.deepcopy() if ok else el for el, ok in zip(self, check))
        return Lattice(elem_generator, elems, iterator=self.attrs_filter, **kwargs)

    # Obsolete methods kept for compatibility
    def radiation_on(self, *args, **kwargs) -> Lattice | None:
        """Obsolete. Turn longitudinal motion on

        The function name is misleading, since the function deals with
        longitudinal motion in general.

        For this reason **the method is obsolete** and replaced by
        :py:meth:`.enable_6d`

        See Also:
            :py:meth:`.enable_6d`
        """
        kwargs.update(zip(("cavity_pass", "dipole_pass", "quadrupole_pass"), args))
        return self._set_6d(True, **kwargs)

    def radiation_off(self, *args, **kwargs) -> Lattice | None:
        """Obsolete. Turn longitudinal motion off

        The function name is misleading, since the function deals with
        longitudinal motion in general.

        For this reason **the method is obsolete** and replaced by
        :py:meth:`.disable_6d`

        See Also:
            :py:meth:`disable_6d`
        """
        kwargs.update(zip(("cavity_pass", "dipole_pass", "quadrupole_pass"), args))
        return self._set_6d(False, **kwargs)

    @property
    def radiation(self) -> bool:
        """Obsolete (see :py:attr:`is_6d` instead)

        :meta private:
        """
        try:
            return self._radiation
        except AttributeError:
            radiate = False
            for elem in self:
                if elem.longt_motion:
                    radiate = True
                    break
            # noinspection PyAttributeOutsideInit
            self._radiation = radiate
            return radiate


def lattice_filter(params, lattice):
    """Copy lattice parameters and run through all lattice elements

    Parameters:
        params:     Dictionary of Lattice attributes
        lattice:    Input ``Lattice``

    Yields:
        lattice ``Elements``
    """
    return lattice.attrs_filter(params, elem_generator, lattice)


# noinspection PyUnusedLocal
def elem_generator(params, elems: Iterable[Element]) -> Iterable[Element]:
    """Run through all elements without any check

    Parameters:
        params:     Dictionary of Lattice attributes
        elems:      Iterable of lattice ``Elements``

    Yields:
        lattice ``Elements``
    """
    return elems


no_filter: Filter = elem_generator  # provided for backward compatibility


def type_filter(params, elems: Iterable[Element]) -> Generator[Element, None, None]:
    """Run through all elements and check element validity.
    Analyse elements for radiation state

    Parameters:
        params:     Dictionary of Lattice attributes
        elems:      Iterable of lattice ``Elements``

    Yields:
        lattice ``Elements``
    """
    radiate = False
    for idx, elem in enumerate(elems):
        if isinstance(elem, Element):
            if elem.longt_motion:
                radiate = True
            yield elem
        else:
            warnings.warn(
                AtWarning(f"item {idx} ({elem}) is not an AT element: " "ignored"),
                stacklevel=2,
            )
    params["_radiation"] = radiate


def params_filter(params, elem_filter: Filter, *args) -> Generator[Element, None, None]:
    """Run through all elements, looking for energy and periodicity.
    Remove the Energy attribute of non-radiating elements

    Parameters:
        params:         Dictionary of Lattice attributes
        elem_filter:    Next ``Elements`` filter
        args:           Arguments forwarded to **elem_filter**

    Yields:
        lattice ``Elements``

    The following keys in ``params`` are set:

    * ``_length``
    * ``periodicity``
    * ``energy`` (optional)
    * ``_frequency`` (optional)
    * ``harmonic_number`` (optional)

    energy is taken from:
        1) The params dictionary
        2) Cavity elements
        3) Any other element

    periodicity is taken from:
        1) The params dictionary
        2) if periodicity <= 0, from the sum of the bending angles of magnets
    """
    el_energies = []
    thetas = []
    cavities = []
    cell_length = 0

    for elem in elem_filter(params, *args):
        if isinstance(elem, elt.RFCavity):
            cavities.append(elem)
        elif hasattr(elem, "Energy"):
            el_energies.append(elem.Energy)
            del elem.Energy
        if isinstance(elem, elt.Dipole):
            thetas.append(elem.BendingAngle)
        cell_length += getattr(elem, "Length", 0.0)
        yield elem

    params["_length"] = cell_length
    cav_energies = [el.Energy for el in cavities if hasattr(el, "Energy")]
    if cavities:
        cavities.sort(key=lambda el: el.Frequency)
        c0 = cavities[0]
        params.setdefault("harmonic_number", getattr(c0, "HarmNumber", math.nan))
        # params['_harmnumber'] = getattr(c0, 'HarmNumber', math.nan)
        params["_frequency"] = getattr(c0, "Frequency", None)

    if "energy" not in params:
        energies = cav_energies or el_energies
        if energies:
            # Guess energy from the Energy attribute of the elements
            energy = params.setdefault("energy", max(energies))
            if min(energies) < max(energies):
                warnings.warn(
                    AtWarning(
                        "Inconsistent energy values, " f'"energy" set to {energy}'
                    ),
                    stacklevel=2,
                )

    if params.setdefault("periodicity", 1) <= 0:
        # Guess periodicity from the bending angles of the lattice
        try:
            nbp = 2.0 * np.pi / sum(thetas)
        except ZeroDivisionError:
            periodicity = 1
            warnings.warn(
                AtWarning('No bending in the cell, set "Periodicity" to 1'),
                stacklevel=2,
            )
        else:
            periodicity = int(round(nbp))
            if abs(periodicity - nbp) > _TWO_PI_ERROR:
                warnings.warn(
                    AtWarning(f"Non-integer number of cells: {nbp} -> {periodicity}"),
                    stacklevel=2,
                )
        params["periodicity"] = periodicity


Lattice.get_uint32_index = get_uint32_index
Lattice.get_bool_index = get_bool_index
Lattice.refcount = _refcount
Lattice.set_shift = set_shift
Lattice.set_tilt = set_tilt
Lattice.set_rotation = set_rotation
Lattice.get_elements = get_elements
Lattice.get_s_pos = get_s_pos
Lattice.select = refpts_iterator
Lattice.get_value_refpts = get_value_refpts
Lattice.set_value_refpts = set_value_refpts
Lattice.get_geometry = get_geometry
