from __future__ import annotations

__all__ = ["RDTObservable"]

from functools import partial

import numpy as np

from .observables import Need, ElementObservable
from ..lattice import Refpts
from ..physics import RDTType

RDT_names = {
    RDTType.FOCUSING: ("h20000", "h00200"),
    RDTType.COUPLING: ("h10010", "h10100"),
    RDTType.CHROMATIC: ("h11001", "h00111", "h20001", "h00201", "h10002"),
    RDTType.GEOMETRIC1: ("h21000", "h30000", "h10110", "h10020", "h10200"),
    RDTType.GEOMETRIC2: (
        "h22000",
        "h11110",
        "h00220",
        "h31000",
        "h40000",
        "h20110",
        "h11200",
        "h20020",
        "h20200",
        "h00310",
        "h00400",
    ),
    RDTType.TUNESHIFT: ("dnux_dJx", "dnux_dJy", "dnuy_dJy"),
}
RDT_code = {nm: code for code, names in RDT_names.items() for nm in names}


def _rdt_access(param, data):
    return getattr(data, param)


class RDTObservable(ElementObservable):
    """Observe a resonance driving term at selected locations.

    Process the local output of :py:func:`.get_rdts`.
    """

    _rdt_type: RDTType

    def __init__(
        self,
        refpts: Refpts,
        param: str,
        name: str | None = None,
        second_order: bool = False,
        **kwargs,
    ):
        # noinspection PyUnresolvedReferences
        r"""Args:
            refpts:         Observation points.
              See ":ref:`Selecting elements in a lattice <refpts>`"
            param:          :ref:`RDT name <rdt_param>`
            name:           Observable name. If :py:obj:`None`, an explicit
              name will be generated
            second_order:   Compute second order terms. Computation is significantly
              longer using this method

        Keyword Args:
            procfun:        Post-processing function. It can be any numpy ufunc or a
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

        The *target*, *weight* and *bounds* inputs must be broadcastable to the
        shape of *value*.

        .. _rdt_param:
        .. rubric:: RDT parameter names

        **h20000**          at.RDTType.FOCUSING
        **h00200**

        **h10010**          at.RDTType.COUPLING
        **h10100**

        **h11001**          at.RDTType.CHROMATIC
        **h00111**
        **h20001**
        **h00201**
        **h10002**

        **h21000**          at.RDTType.GEOMETRIC1
        **h30000**
        **h10110**
        **h10020**
        **h10200**

        **h22000**          at.RDTType.GEOMETRIC2
        **h11110**
        **h00220**
        **h31000**
        **h40000**
        **h20110**
        **h11200**
        **h20020**
        **h20200**
        **h00310**
        **h00400**

        **dnux_dJx**        at.RDTType.DETUNING
        **dnux_dJy**
        **dnuy_dJy**

        Examples:

            >>> obs = RDTObservable(at.Monitor, "h20000")

            Observe "h20000" at all :py:class:`.Monitor` locations
        """

        needs = {Need.RDT}
        if second_order:
            needs.add(Need.RDT_2ND_ORDER)
        name = self._set_name(name, param, None)
        fun = partial(_rdt_access, param)

        super().__init__(fun, refpts, needs=needs, name=name, **kwargs)
        self._rdt_type = RDT_code[param]
