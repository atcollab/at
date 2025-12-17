__all__ = ["LocalEmittanceObservable"]

from collections.abc import Callable
from functools import partial
from typing import ClassVar

import numpy as np

from ..lattice import AxisDef, plane_, axis_
from ..lattice import Refpts
from .observables import ElementObservable, Need
# noinspection PyProtectedMember
from .observables import _all_rows, _record_access, _fun_access, _subscript


def _sigma4(_, index, data, **kwargs):
    """Betatron beam size: square root of diagonal terms of r44."""
    idx = index[1]
    return np.sqrt(data.r44[_all_rows((idx, idx))])


def _sigma6(_, index, data, **kwargs):
    """Beam size: square root of diagonal terms of r66."""
    idx = index[1]
    return np.sqrt(data.r66[_all_rows((idx, idx))])


_opdata = {"sigma4": _sigma4, "sigma6": _sigma6}


class LocalEmittanceObservable(ElementObservable):
    """Observe emittance-related parameters."""

    # Class attributes
    _pinfo: ClassVar[dict] = {
        "r66": (r"$R66_{{{plane}}}$", "R", _subscript),
        "r44": (r"$R44_{{{plane}}}$", "R", _subscript),
        "emitXY": (
            r"$\epsilon_{{{plane}}}$",
            "Emittance [m]",
            partial(plane_, key="label"),
        ),
        "emitXYZ": (
            r"$\epsilon_{{{plane}}}$",
            "Emittance [m]",
            partial(plane_, key="label"),
        ),
        "m66": (r"$\mathrm{{T}}_{{{plane}}}$", "T [m]", _subscript),
        "orbit6": (r"${{{plane}}}_{{co}}$", "closed orbit", partial(axis_, key="code")),
        "sigma4": (r"$\sigma_{{{plane}}}$", "sigma", partial(axis_, key="label")),
        "sigma6": (r"$\sigma_{{{plane}}}$", "sigma", partial(axis_, key="label")),
    }

    def __init__(
        self,
        refpts: Refpts,
        param: str | Callable,
        *,
        plane: AxisDef = Ellipsis,
        axis: AxisDef = Ellipsis,
        name: str | None = None,
        summary: bool = False,
        label: str | None = None,
        **eval_kw,
    ):
        # noinspection PyUnresolvedReferences
        r"""Args:
            refpts:     Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            param:      :ref:`Emittance parameter name <localemit_param>`
              or :ref:`user-defined evaluation function <localemit_eval>`
            axis:       Index in the parameter array. If :py:obj:`Ellipsis`, the
              whole array is specified
            plane:      Index in the parameter array for *emitXY* and *emitXYZ*
              parameters. If :py:obj:`Ellipsis`, the whole array is specified
            name:       Observable name. If :py:obj:`None`, an explicit
              name will be generated.

        Keyword Args:
            summary:        Set to :py:obj:`True` if the user-defined
             evaluation function returns a single item (see below)
            postfun:        Post-processing function. It can be any numpy ufunc or a
              function name in {"real", "imag", "abs", "angle", "log", "exp", "sqrt"}.
            statfun:        Statistics post-processing function. it can be a numpy
              function or a function name in {"mean", "std", "var", "min", "max"}.
              Example: :pycode:`statfun=numpy.mean`.
            target:         Target value for a constraint. If :py:obj:`None`
              (default), the residual will always be zero.
            weight:         Weight factor: the residual is
              :pycode:`((value-target)/weight)**2`
            bounds:         Tuple of lower and upper bounds. The parameter
              is constrained in the interval
              [*target*\ +\ *low_bound* *target*\ +\ *up_bound*]

        .. rubric:: Evaluation keywords

        These values must be provided to the :py:meth:`~.ObservableList.evaluate`
        method. Default values may be given at instantiation.

        * **ring**:     Lattice description,
        * **dp**:     Momentum deviation. Defaults to :py:obj:`None`,
        * **dct**:    Path lengthening. Defaults to :py:obj:`None`,
        * **df**:     Deviation from the nominal RF frequency.
          Defaults to :py:obj:`None`,
        * **orbit**:   Initial orbit. Avoids looking for the closed orbit if it is
          already known,

        .. rubric:: Shape of the value

        If the requested attribute has shape :pycode:`shp`, then
        *value* has shape :pycode:`(nrefs,) + shp` with :pycode:`nrefs` number of
        reference points: a  :pycode:`nrefs` dimension is prepended to the shape of
        the attribute, **even if nrefs == 1**. The *target*, *weight* and *bounds*
        inputs must be broadcastable to the shape of *value*. For instance, a *target*
        with shape :pycode:`shp` will automatically broadcast and apply  to all
        reference points.

        .. _localemit_param:
        .. rubric:: Emittance parameter names

        In addition to :py:func:`.ohmi_envelope` parameter names,
        LocalEmittanceObservable adds the *sigma* parameter:

        ===========    ===================================================
        **r66**        (6, 6) equilibrium envelope matrix R
        **r44**        (4, 4) betatron emittance matrix (dp = 0)
        **m66**        (6, 6) transfer matrix from the start of the ring
        **orbit6**     (6,) closed orbit
        **emitXY**     (2,) betatron emittance projected on xxp and yyp
        **emitXYZ**    (3,) 6x6 emittance projected on xxp, yyp, ldp
        **sigma**      (6,) standard deviation of the beam (square root of the diagonal
                            terms of *r66*)
        ===========    ===================================================

        .. _localemit_eval:
        .. rubric:: User-defined evaluation function

        The observable value is computed as:

        :pycode:`value = fun(emit, **kwargs)[plane]`

        - *emit* is the output of :py:func:`.ohmi_envelope`, evaluated at the
          *refpts* of the observable,
        - *kwargs* are the keyword arguments provided to the observable constructor,
          to the constructor of the enclosing :py:class:`.ObservableList` and to the
          :py:meth:`~.ObservableList.evaluate` method,
        - *value* is the value of the Observable and must have one line per
          refpoint. Alternatively, it may be a single line, but then the
          *summary* keyword must be set to :py:obj:`True`,
        - the *plane* or *axis* keyword then selects the desired values in the function
          output.

        Examples:
            >>> obs = LocalEmittanceObservable(at.Monitor, "sigma")

            Observe the beam standard deviation in all 6 axes at Monitor locations

            >>> obs = LocalEmittanceObservable(
            ...     at.Quadrupole, "sigma", axis="x", statfun=np.max
            ... )

            Observe the maximum horizontal beam size in Quadrupoles
        """
        needs = {Need.LOCALEMIT}
        descr = plane_(plane) if param in {"emitXY", "emitXYZ"} else axis_(axis)
        name = self._set_name(name, param, descr["code"])
        index = _all_rows(descr["index"])
        if callable(param):
            if summary:
                fun = partial(_fun_access, param, descr["index"])
            else:
                fun = partial(_fun_access, param, index)
        else:
            fun = partial(_opdata.get(param, _record_access), param, index)

        super().__init__(
            fun,
            refpts,
            needs=needs,
            name=name,
            summary=summary,
            label=self._pl_lab(param, descr["index"]),
            axis_label=self._ax_lab(param, descr["index"]),
            **eval_kw,
        )
        if label:
            self.label = label
