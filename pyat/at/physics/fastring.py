"""
Functions relating to fast_ring.
"""

from __future__ import annotations

__all__ = ["fast_ring", "simple_ring"]

from collections.abc import Sequence
import numpy as np
from ..lattice import Lattice, Refpts, AtError, Particle
from ..lattice import Drift, RFCavity, Element, Marker, SimpleQuantDiff
from ..lattice import DeltaQ, SimpleRadiation, M66
from ..physics import gen_m66_elem, gen_detuning_elem, gen_quantdiff_elem
from ..constants import clight


def _replace_cav(cav_element: Element) -> Element:
    name = ("DR_" + cav_element.FamName,)
    elem = Drift(name, cav_element.Length) if cav_element.Length > 0 else Marker(name)
    return elem


def _merge_cavs(all_cavs: Sequence) -> Sequence[RFCavity]:
    m_cavs = []
    freqs = [e.Frequency for e in all_cavs]
    for fr in np.atleast_1d(np.unique(freqs)):
        cavf = [cav for cav in all_cavs if cav.Frequency == fr]
        vol = np.sum([c.Voltage for c in cavf])
        cavl = RFCavity(
            cavf[0].FamName,
            0,
            vol,
            fr,
            cavf[0].HarmNumber,
            cavf[0].Energy,
            TimeLag=cavf[0].TimeLag,
        )
        m_cavs.append(cavl)
    return m_cavs


def _rearrange(ring) -> tuple:
    cav_idx = ring.get_uint32_index(RFCavity)
    mcavs = _merge_cavs(ring[cav_idx])
    ring[cav_idx] = [_replace_cav(ring[i]) for i in cav_idx]
    ring += mcavs
    ring += [Marker("xsplit")]


def _split_ring(
    ring: Lattice, split_inds: Refpts | None = None, keep_cavities: bool = False
) -> Sequence:
    inds = ring.get_bool_index(split_inds, endpoint=True)
    inds[[0, -1]] = True
    if keep_cavities:
        cavsb = ring.get_uint32_index(RFCavity)
        cavse = cavsb + 1
        inds = inds | ring.get_bool_index(cavsb) | ring.get_bool_index(cavse)
    inds = ring.get_uint32_index(inds)
    split_ring = []
    for b, e in zip(inds[:-1], inds[1:], strict=True):
        split_ring += [Marker("xsplit")]
        split_ring += ring[int(b) : int(e)]
    split_ring += [Marker("xsplit")]
    split_ring = Lattice(split_ring, **vars(ring))
    if keep_cavities is False:
        _rearrange(split_ring)
    return split_ring, split_ring.get_uint32_index("xsplit")


def fast_ring(
    ring,
    split_inds: Refpts | None = None,
    qpx: Sequence[float] | None = None,
    qpy: Sequence[float] | None = None,
    detuning_coeff: Sequence[float] | None = None,
    alphac: Sequence[float] | None = None,
    emitx: float | None = None,
    emity: float | None = None,
    espread: float | None = None,
    keep_cavities: bool = False,
    dp: float | None = None,
) -> Lattice:
    """
    A fast ring consisting in the following elements.

    * a RF cavity per distinct frequency,
    * a 6x6 linear transfer map,
    * a detuning and chromaticity element,
    * a quantum diffusion element active when radiations are enabled.

    A new lattices is returned with the same attributes (energy, particle,
    circumference, periodicity,…) as the initial one. The radiation state is
    also concerved. Radiations can be turned on or of using ring.disable_6d()
    or ring.enable_6d(). *Rad attribute will be used when radiations are turned
    on.

    It is possible to split the original ring in multiple "fastrings"
    using the ``split_inds`` argument.

    By default all the cavities are merged ans placed at the end of the ring this
    allows to minimize the number of segments and therefore lattice elements.
    Howevever, the longitudinal dynamics is modified by this action. In order to
    preserve the original RF configuration the ``keep_cavities`` argument should
    be set to true. In this configuration the tracking is slower as more segment are
    introduced.

    Parameters:
        ring:       Lattice description
        split_inds: List of indexes where to split the ring
        qpx:            horizontal chromatic detuning coefficients. Default None.
          If specified qpy should also be provided, if not first order term computed
          automatically from ring. at.physics.chromaticity can be used to compute
          higher order chromatic coefficients
        qpy:            vertical chromatic detuning coefficients. Default None.
          If specified qpy should also be provided, if not first order term computed
          automatically from ring. at.physics.chromaticity can be used to compute
          higher order chromatic coefficients
        detuning_coeff: First order amplitude detuning coefficients
          [dQx/dJx, dQx/dJy, dQy/dJy]. Default None: coefficient computer from ring.
          at.physics.detuning can be used to compute detuning coefficients
        alphac: Higher order momentum compaction factor. Default None. Can be
          Can be claculated with physics.get_mcf()
        emitx: Equilibrium horizontal emittance. Default is the one given by ring
        emity: Equilibrium vertical emittance. Default is the one given by ring
        espread: Equilibrium energy spread. Default is the one given by ring
        keep_cavities: The ring will be sliced at each cavities. This allows to
          match longitudinal dynamics of the original ring bu t slows down the
          tracking as more segments are generated

    Returns:
        fastring (Lattice):    Fast ring lattice object
    """
    split_ring, split_inds = _split_ring(ring, split_inds, keep_cavities)
    if dp is not None:
        split_ring.set_rf_frequency(dp=dp)
    _, o4 = split_ring.disable_6d(copy=True).find_orbit(refpts=split_inds, dp=dp)
    _, o6 = split_ring.enable_6d(copy=True).find_orbit(refpts=split_inds, dp=dp)

    fastring = []
    for o4b, o4e, o6b, o6e, sib, sie in zip(
        o4[:-1], o4[1:], o6[:-1], o6[1:], split_inds[:-1], split_inds[1:], strict=True
    ):
        r = split_ring[sib + 1 : sie]
        iscav = np.array([isinstance(elem, RFCavity) for elem in r], dtype=bool)
        if np.all(iscav):
            fastring = [*fastring, *list(r)]
        elif np.all(~iscav):
            lin_elem = gen_m66_elem(
                r.disable_6d(copy=True), o4b, o4e, r.enable_6d(copy=True), o6b, o6e
            )
            fastring = [*fastring, lin_elem]
        else:
            msg = (
                "One fast ring segment contains a combination of magnets and cavity."
                "They need to be separated"
            )
            raise AtError(msg)
    detuning_elem = gen_detuning_elem(
        split_ring,
        qpx=qpx,
        qpy=qpy,
        detuning_coeff=detuning_coeff,
        alphac=alphac,
        orbit=o4[-1],
        orbit6=o6[-1],
    )
    if emity is None and emitx is None and espread is None:
        qd_elem = gen_quantdiff_elem(split_ring.enable_6d(copy=True), orbit=o6[-1])
    else:
        params = split_ring.enable_6d(copy=True).envelope_parameters()
        taux, tauy, tauz = params.Tau
        if emitx is None:
            emitx = params.emittances[0]
        if emity is None:
            emity = params.emittances[1]
        if espread is None:
            espread = params.sigma_e
        qd_elem = SimpleQuantDiff(
            "Diffusion",
            betax=detuning_elem.BetaRad[0],
            betay=detuning_elem.BetaRad[1],
            emitx=emitx,
            emity=emity,
            espread=espread,
            taux=taux * clight / ring.circumference,
            tauy=tauy * clight / ring.circumference,
            tauz=tauz * clight / ring.circumference,
        )
    fastring.append(detuning_elem)
    fastring.append(qd_elem)
    fastring = Lattice(fastring, **vars(ring))
    fastring.periodicity = 1
    if ring.radiation:
        fastring.enable_6d()
    else:
        fastring.disable_6d()
    return fastring


def simple_ring(
    energy: float,
    circumference: float,
    harmonic_number: float | Sequence[float],
    Qx: float,
    Qy: float,
    Vrf: float | Sequence[float],
    alpha: float,
    betax: float = 1.0,
    betay: float = 1.0,
    alphax: float = 0.0,
    alphay: float = 0.0,
    dispx: float = 0.0,
    dispxp: float = 0.0,
    dispy: float = 0.0,
    dispyp: float = 0.0,
    Qpx: float | Sequence[float] = 0.0,
    Qpy: float | Sequence[float] = 0.0,
    A1: float = 0.0,
    A2: float = 0.0,
    A3: float = 0.0,
    emitx: float = 0.0,
    emity: float = 0.0,
    espread: float = 0.0,
    taux: float = 0.0,
    tauy: float = 0.0,
    tauz: float = 0.0,
    U0: float = 0.0,
    name: str = "",
    particle: str | Particle = "relativistic",
    TimeLag: float | Sequence[float] = 0.0,
) -> Lattice:
    """Generates a "simple ring" based on a given dictionary of global parameters.

    A simple ring consists of:

    * an RF cavity,
    * a simple radiation damping element
    * a simplified quantum diffusion element which contains equilibrium emittance
    * a detuning and chromaticity element,
    * a 6x6 linear transfer map with no radiation damping

    Parameters:
        energy: [eV]
        circumference: [m]
        harmonic_number: can be scalar or sequence of scalars. The RF
          frequency is derived from this and the ring circumference
        Qx: horizontal tune
        Qy: vertical tune
        Vrf: RF Voltage set point [V] - can be scalar or sequence of scalars
        alpha: momentum compaction factor
        betax: horizontal beta function [m], Default=1
        betay: vertical beta function [m], Default=1
        alphax: horizontal alpha function, Default=0
        alphay: vertical alpha function, Default=0
        dispx: horizontal dispersion [m], Default=0
        dispxp: horizontal dispersion prime, Default=0
        dispy: vertical dispersion [m], Default=0
        dispyp: vertical dispersion prime, Default=0
        Qpx: If single value, it is horizontal linear chromaticity
          If an array is given it corresponds to a list of horizontal
          non linear chromaticities [Q',Q'',Q''',...]. This is expanded
          following Q'/1! * (dp/p) + Q''/2! *(dp/p)^2 etc. Default=0.0
        Qpy: If single value, it is vertical linear chromaticity
          If an array is given it corresponds to a list of horizontal
          non linear chromaticities [Q',Q'',Q''',...]. This is expanded
          following Q'/1! * (dp/p) + Q''/2! *(dp/p)^2 etc. Default=0.0
        A1: horizontal amplitude detuning coefficient, Default=0
        A2: cross term for amplitude detuning coefficient, Default=0
        A3: vertical amplitude detuning coefficient, Default=0
        emitx: horizontal equilibrium emittance [m.rad], Default=0
          ignored if emitx=0
        emity: vertical equilibrium emittance [m.rad], Default=0
          ignored if emity=0
        espread: equilibrium momentum spread, Default=0
          ignored if espread=0
        taux: horizontal radiation damping time [turns], Default=0
          ignored if taux=0
        tauy: vertical radiation damping time [turns], Default=0
          ignored if tauy=0
        tauz: longitudinal radiation damping time [turns], Default=0
          ignored if tauz=0
        U0: energy loss [eV] (positive number), Default=0
        name: Name of the lattice
        particle: circulating particle. May be
          'relativistic', 'electron', 'positron', 'proton'
          or a Particle object
        TimeLag: Set the timelag of the cavities, Default=0. Can be scalar
          or sequence of scalars (as with harmonic_number and Vrf).

    If the given emitx, emity or espread is 0, then no equlibrium emittance
    is applied in this plane.
    If the given tau is 0, then no radiation damping is applied for this plane.

    Returns:
        ring:    Simple ring
    """
    try:
        rfp = np.broadcast(Vrf, harmonic_number, TimeLag)
    except ValueError as exc:
        msg = "Vrf, harmonic_number and TimeLag must be broadcastable"
        raise AtError(msg) from exc

    # revolution frequency
    f0 = clight / circumference

    all_cavities = [
        RFCavity(f"RFC{i+1}", 0.0, v, h * f0, h, energy, TimeLag=t)
        for i, (v, h, t) in enumerate(rfp)
    ]

    # Now we will use the optics parameters to compute the uncoupled M66 matrix

    s_dphi_x = np.sin(2 * np.pi * Qx)
    c_dphi_x = np.cos(2 * np.pi * Qx)
    s_dphi_y = np.sin(2 * np.pi * Qy)
    c_dphi_y = np.cos(2 * np.pi * Qy)

    M00 = c_dphi_x + alphax * s_dphi_x
    M01 = betax * s_dphi_x
    M10 = -(1.0 + alphax**2) / betax * s_dphi_x
    M11 = c_dphi_x - alphax * s_dphi_x

    M04 = (1 - M00) * dispx - M01 * dispxp
    M14 = -M10 * dispx + (1 - M11) * dispxp

    M22 = c_dphi_y + alphay * s_dphi_y
    M23 = betay * s_dphi_y
    M32 = -(1.0 + alphay**2) / betay * s_dphi_y
    M33 = c_dphi_y - alphay * s_dphi_y

    M24 = (1 - M22) * dispy - M23 * dispyp
    M34 = -M32 * dispy + (1 - M33) * dispyp

    M44 = 1.0
    M45 = 0.0
    M54 = alpha * circumference
    M55 = 1

    Mat66 = np.array(
        [
            [M00, M01, 0.0, 0.0, M04, 0.0],
            [M10, M11, 0.0, 0.0, M14, 0.0],
            [0.0, 0.0, M22, M23, M24, 0.0],
            [0.0, 0.0, M32, M33, M34, 0.0],
            [0.0, 0.0, 0.0, 0.0, M44, M45],
            [0.0, 0.0, 0.0, 0.0, M54, M55],
        ],
        order="F",
    )

    # generate the linear tracking element, we set a length
    # which is needed to give the lattice object the correct length
    # (although it is not used for anything else)
    lin_elem = M66("Linear", m66=Mat66, Length=circumference)

    # Generate the simple radiation element
    simplerad = SimpleRadiation(
        "SR",
        taux=taux,
        tauy=tauy,
        tauz=tauz,
        U0=U0,
        dispx=dispx,
        dispy=dispy,
        dispxp=dispxp,
        dispyp=dispyp,
    )

    # Generate the simple quantum diffusion element
    quantdiff = SimpleQuantDiff(
        "SQD",
        betax=betax,
        betay=betay,
        emitx=emitx,
        emity=emity,
        espread=espread,
        taux=taux,
        tauy=tauy,
        tauz=tauz,
    )

    chromx_arr = np.ravel(Qpx)
    chromy_arr = np.ravel(Qpy)
    chrom_maxorder = max(chromx_arr.size, chromy_arr.size)

    chromx_arr = np.pad(chromx_arr, (0, chrom_maxorder - len(chromx_arr)))
    chromy_arr = np.pad(chromy_arr, (0, chrom_maxorder - len(chromy_arr)))

    # Generate the detuning element
    nonlin_elem = DeltaQ(
        "NonLinear",
        alpha=[alphax, alphay],
        beta=[betax, betay],
        qpx=chromx_arr,
        qpy=chromy_arr,
        detuning_coefficients=[A1, A2, A3],
        Chrom_MaxOrder=chrom_maxorder,
    )

    # Assemble all elements into the lattice object
    ring = Lattice(
        [*all_cavities, simplerad, quantdiff, lin_elem, nonlin_elem],
        name=name,
        energy=energy,
        particle=particle,
        periodicity=1,
    )

    return ring
