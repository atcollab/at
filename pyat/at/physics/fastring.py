"""
Functions relating to fast_ring
"""

from __future__ import annotations

__all__ = ["fast_ring", "simple_ring"]

import copy
from collections.abc import Sequence
from functools import reduce

import numpy as np

from ..constants import clight
from ..lattice import Lattice, Particle, Refpts
from ..lattice import RFCavity, Element, Marker, get_cells, checkname
from ..lattice import get_elements, M66, SimpleQuantDiff, AtError, SimpleRadiation
from ..physics import gen_m66_elem, gen_detuning_elem, gen_quantdiff_elem


def _rearrange(ring: Lattice, split_inds: Refpts = None):
    inds = ring.get_bool_index(split_inds, endpoint=True)
    inds[0] = True
    inds[-1] = True
    inds = ring.get_uint32_index(inds)
    all_rings = [ring[int(b) : int(e)] for b, e in zip(inds[:-1], inds[1:])]

    ringm = []
    for ring_slice in all_rings:
        ring_slice.insert(0, Marker("xbeg"))
        ring_slice.append(Marker("xend"))
        cavs = [e for e in ring_slice if isinstance(e, RFCavity)]
        newpass = ["IdentityPass" if c.Length == 0 else "DriftPass" for c in cavs]
        for c, pm in zip(cavs, newpass):
            c.PassMethod = pm
        uni_freq = np.unique([e.Frequency for e in cavs])
        for fr in np.atleast_1d(uni_freq):
            cavf = [c for c in cavs if c.Frequency == fr]
            vol = reduce(lambda x, y: x + y, (c.Voltage for c in cavf))
            cavl = RFCavity("CAVL", 0, vol, fr, cavf[0].HarmNumber, cavf[0].Energy)
            cavl.TimeLag = cavf[0].TimeLag
            ring_slice.append(cavl)
        ringm = ringm + ring_slice
    return all_rings, Lattice(ringm, energy=ring.energy, periodicity=1)


def _fring(ring, split_inds: Refpts = None, detuning_elem=None):
    all_rings, merged_ring = _rearrange(ring, split_inds=split_inds)
    ibegs = get_cells(merged_ring, checkname("xbeg"))
    iends = get_cells(merged_ring, checkname("xend"))
    _, orbit = merged_ring.find_orbit(refpts=ibegs | iends)
    if detuning_elem is None:
        detuning_elem = gen_detuning_elem(merged_ring, orbit[-1])
    else:
        detuning_elem.T1 = -orbit[-1]
        detuning_elem.T2 = orbit[-1]

    fastring = []
    for counter, r in enumerate(all_rings):
        cavs = [e for e in r if e.PassMethod.endswith("CavityPass")]
        [r.remove(c) for c in cavs]
        lin_elem = gen_m66_elem(r, orbit[2 * counter], orbit[2 * counter + 1])
        lin_elem.FamName = lin_elem.FamName + "_" + str(counter)
        [fastring.append(c) for c in cavs]
        fastring.append(lin_elem)
    fastring.append(detuning_elem)
    try:
        qd_elem = gen_quantdiff_elem(merged_ring)
        fastring.append(qd_elem)
    except ValueError:  # No synchrotron radiation => no diffusion element
        pass
    fastring = Lattice(fastring, **vars(ring))
    return fastring


def fast_ring(ring: Lattice, split_inds: Refpts = None) -> tuple[Lattice, Lattice]:
    """Generates a "fast ring"

    A fast ring consisting in:

    * a RF cavity per distinct frequency,
    * a 6x6 linear transfer map,
    * a detuning and chromaticity element,
    * a quantum diffusion element (for radiation ring).


    2 new lattices are returned, one with radiation and one without. These lattices
    keep the same attributes (energy, particle, circumference, periodicity,…)
    as the initial one.

    It is possible to split the original ring in multiple "fastrings"
    using the ``split_inds`` argument.

    Parameters:
        ring:       Lattice description
        split_inds: List of indexes where to split the ring

    Returns:
        fring (Lattice):    Fast ring without radiation
        fringrad (Lattice): Fast ring with radiation
    """
    ringi = ring.deepcopy()
    fastringnorad = _fring(ringi.radiation_off(copy=True), split_inds=split_inds)
    detuning_elem = copy.deepcopy(get_elements(fastringnorad, "NonLinear")[0])
    fastringrad = _fring(
        ringi.radiation_on(copy=True),
        split_inds=split_inds,
        detuning_elem=detuning_elem,
    )
    return fastringnorad, fastringrad


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
    TimeLag: float | Sequence[float] = 0.0
) -> Lattice:
    """Generates a "simple ring" based on a given dictionary
       of global parameters

    A simple ring consists of:

    * an RF cavity,
    * a 6x6 linear transfer map with no radiation damping,
    * a detuning and chromaticity element,
    * a simple radiation damping element
    * a simplified quantum diffusion element which contains equilibrium emittance

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
        raise AtError("Vrf, harmonic_number and TimeLag must be broadcastable") from exc

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

    chromx_arr = np.pad(chromx_arr, (0, chrom_maxorder-len(chromx_arr)))
    chromy_arr = np.pad(chromy_arr, (0, chrom_maxorder-len(chromy_arr)))

    # Generate the detuning element
    nonlin_elem = Element(
        "NonLinear",
        PassMethod="DeltaQPass",
        Betax=betax,
        Betay=betay,
        Alphax=alphax,
        Alphay=alphay,
        chromx_arr=chromx_arr,
        chromy_arr=chromy_arr,
        A1=A1,
        A2=A2,
        A3=A3,
        chrom_maxorder=chrom_maxorder,
    )

    # Assemble all elements into the lattice object
    ring = Lattice(
        all_cavities + [lin_elem, nonlin_elem, simplerad, quantdiff],
        name=name,
        energy=energy,
        particle=particle,
        periodicity=1,
    )

    return ring
