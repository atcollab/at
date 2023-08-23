"""
Functions relating to fast_ring
"""
from functools import reduce
import numpy
from typing import Tuple
from at.lattice import RFCavity, Element, Marker, Lattice, get_cells, checkname
from at.lattice import get_elements, M66, EnergyLoss, SimpleQuantDiff
from at.physics import gen_m66_elem, gen_detuning_elem, gen_quantdiff_elem
from at.constants import clight, e_mass
import copy


__all__ = ['fast_ring', 'gen_simple_ring']


def _rearrange(ring: Lattice, split_inds=[]):
    inds = numpy.append(split_inds, [0, len(ring)+1])
    inds = numpy.unique(inds)
    all_rings = [ring[int(b):int(e)] for b, e in zip(inds[:-1], inds[1:])]

    ringm = []
    for ring_slice in all_rings:
        ring_slice.insert(0, Marker('xbeg'))
        ring_slice.append(Marker('xend'))
        cavs = [e for e in ring_slice if isinstance(e, RFCavity)]
        newpass = ['IdentityPass' if c.Length == 0
                   else 'DriftPass' for c in cavs]
        for c, pm in zip(cavs, newpass):
            c.PassMethod = pm
        uni_freq = numpy.unique([e.Frequency for e in cavs])
        for fr in numpy.atleast_1d(uni_freq):
            cavf = [c for c in cavs if c.Frequency == fr]
            vol = reduce(lambda x, y: x+y, (c.Voltage for c in cavf))
            cavl = RFCavity('CAVL', 0, vol, fr,
                            cavf[0].HarmNumber, cavf[0].Energy)
            cavl.TimeLag = cavf[0].TimeLag
            ring_slice.append(cavl)
        ringm = ringm + ring_slice
    return all_rings, Lattice(ringm, energy=ring.energy)


def _fring(ring, split_inds=[], detuning_elem=None):
    all_rings, merged_ring = _rearrange(ring, split_inds=split_inds)
    ibegs = get_cells(merged_ring, checkname('xbeg'))
    iends = get_cells(merged_ring, checkname('xend'))
    _, orbit = merged_ring.find_orbit(refpts=ibegs | iends)
    if detuning_elem is None:
        detuning_elem = gen_detuning_elem(merged_ring, orbit[-1])
    else:
        detuning_elem.T1 = -orbit[-1]
        detuning_elem.T2 = orbit[-1]

    fastring = []
    for counter, r in enumerate(all_rings):
        cavs = [e for e in r if e.PassMethod.endswith('CavityPass')]
        [r.remove(c) for c in cavs]
        lin_elem = gen_m66_elem(r, orbit[2*counter],
                                orbit[2*counter+1])
        lin_elem.FamName = lin_elem.FamName + '_' + str(counter)
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


def fast_ring(ring: Lattice, split_inds=[]) -> Tuple[Lattice, Lattice]:
    """Generates a "fast ring"

    A fast ring consisting in:

    * a RF cavity per distinct frequency,
    * a 6x6 linear transfer map,
    * a detuning and chromaticity element,
    * a quantum diffusion element (for radiation ring).

    2 new lattices are returned, one with radiation and one without
    It is possible to split the original ring in multiple "fastrings"
    using the ``split_inds`` argument

    Parameters:
        ring:       Lattice description
        split_inds: List of indexes where to split the ring

    Returns:
        fring (Lattice):    Fast ring without radiation
        fringrad (Lattice): Fast ring with radiation
    """
    ringi = ring.deepcopy()
    fastringnorad = _fring(ringi.radiation_off(copy=True),
                           split_inds=split_inds)
    detuning_elem = copy.deepcopy(get_elements(fastringnorad,
                                               'NonLinear')[0])
    fastringrad = _fring(ringi.radiation_on(copy=True),
                         split_inds=split_inds,
                         detuning_elem=detuning_elem)
    return fastringnorad, fastringrad


def gen_simple_ring(ring_dictionary):
    """Generates a "simple ring" based on a given dictionary
       of global parameters

    A simple ring consists of:

    * an RF cavity,
    * a 6x6 linear transfer map,
    * a detuning and chromaticity element,
    * an energy loss element
    * a simplified quantum diffusion element

    Parameters:
        ring_dictionary:       Lattice description
        The ring_dictionary must contain the following parameters,
            * energy [eV]
            * circumference [m]
            * harmonic_number
            * alpha_x, alpha_y
            * beta_x, beta_y [m]
            * Qx, Qy - full or fractional tunes
            * alpha (momentum compaction factor)
            * U0 - energy loss [eV] (positive number)
            * Vrf - RF Voltage set point [V]
            * Qpx, Qpy - linear chromaticities
            * A1, A2, A3 - amplitude detuning coefficients
            * emit_x, emit_y, sigma_dp - equilibrium values [m.rad, m.rad, -]

    Returns:
        ring (Lattice):    Simple ring
    """

    # parse everything first

    circumference = ring_dictionary['circumference']
    harmonic_number = ring_dictionary['harmonic_number']

    energy = ring_dictionary['energy']

    alpha_x = ring_dictionary['alpha_x']
    alpha_y = ring_dictionary['alpha_y']

    beta_x = ring_dictionary['beta_x']
    beta_y = ring_dictionary['beta_y']

    Qx = ring_dictionary['Qx']
    Qy = ring_dictionary['Qy']

    alpha = ring_dictionary['alpha']
    U0 = ring_dictionary['U0']
    Vrf = ring_dictionary['Vrf']

    Qpx = ring_dictionary['Qpx']
    Qpy = ring_dictionary['Qpy']

    A1 = ring_dictionary['A1']
    A2 = ring_dictionary['A2']
    A3 = ring_dictionary['A3']

    emit_x = ring_dictionary['emit_x']
    emit_y = ring_dictionary['emit_y']
    sigma_dp = ring_dictionary['sigma_dp']

    tau = ring_dictionary['tau']
    
    # compute rf frequency
    frf = harmonic_number * clight / circumference

    # compute slip factor
    gamma = energy / e_mass
    eta = alpha-1/gamma**2

    # compute the synchronous phase and the TimeLag
    phi_s = numpy.arcsin(U0/Vrf)
    TimeLag = clight * phi_s / (2 * numpy.pi * frf)
    # generate rf cavity element
    rfcav = RFCavity('RFC', 0, Vrf, frf, harmonic_number, energy,
                     TimeLag=TimeLag)

    # Now we will use the optics parameters to compute the uncoupled M66 matrix
    def sincos(x):
        return numpy.sin(x), numpy.cos(x)

    s_dphi_x, c_dphi_x = sincos(2*numpy.pi*Qx)
    s_dphi_y, c_dphi_y = sincos(2*numpy.pi*Qy)

    M00 = c_dphi_x + alpha_x * s_dphi_x
    M01 = beta_x * s_dphi_x
    M10 = -(1. + alpha_x**2) / beta_x * s_dphi_x
    M11 = c_dphi_x - alpha_x * s_dphi_x

    M22 = c_dphi_y + alpha_y * s_dphi_y
    M23 = beta_y * s_dphi_y
    M32 = -(1. + alpha_y**2) / beta_y * s_dphi_y
    M33 = c_dphi_y - alpha_y * s_dphi_y

    M44 = 1.
    M45 = 0.
    M54 = eta*circumference
    M55 = 1

    Mat66 = numpy.array([[M00, M01, 0., 0., 0., 0.],
                         [M10, M11, 0., 0., 0., 0.],
                         [0., 0., M22, M23, 0., 0.],
                         [0., 0., M32, M33, 0., 0.],
                         [0., 0., 0., 0., M44, M45],
                         [0., 0., 0., 0., M54, M55]], order='F')

    # generate the linear tracking element, we set a length
    # which is needed to give the lattice object the correct length
    # (although it is not used)
    lin_elem = M66('Linear', m66=Mat66, Length=circumference)

    # Generate the simple quantum diffusion element
    quantdiff = SimpleQuantDiff('SQD', emit_x, emit_y, sigma_dp,
                                tau[0], tau[1], tau[2],
                                beta_x, beta_y, U0)

    # Generate the detuning element
    nonlin_elem = Element('NonLinear', PassMethod='DeltaQPass',
                          Betax=beta_x, Betay=beta_y,
                          Alphax=alpha_x, Alphay=alpha_y,
                          Qpx=Qpx, Qpy=Qpy,
                          A1=A1, A2=A2, A3=A3)

    # Assemble all elements into the lattice object
    ring = Lattice([rfcav, lin_elem, nonlin_elem, quantdiff],
                   energy=energy)

    return ring
