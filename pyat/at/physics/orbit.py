"""
Closed orbit related functions
"""
import numpy
from at.lattice.constants import clight
from at.lattice import AtWarning, check_radiation, DConstant
from at.lattice import Lattice, get_s_pos, uint32_refpts
from at.tracking import lattice_pass
from at.physics import ELossMethod, get_timelag_fromU0
import warnings

__all__ = ['find_orbit4', 'find_sync_orbit', 'find_orbit6', 'find_orbit']


@check_radiation(False)
def _orbit_dp(ring, dp=None, guess=None, **kwargs):
    """Solver for fixed energy deviation"""
    # We seek
    #  - f(x) = x
    #  - g(x) = f(x) - x = 0
    #  - g'(x) = f'(x) - 1
    # Use a Newton-Raphson-type algorithm:
    #  - r_n+1 = r_n - g(r_n) / g'(r_n)
    #  - r_n+1 = r_n - (f(r_n) - r_n) / (f'(r_n) - 1)
    #
    # (f(r_n) - r_n) / (f'(r_n) - 1) can be seen as x = b/a where we use least
    #     squares fitting to determine x when ax = b
    # f(r_n) - r_n is denoted b
    # f'(r_n) is the 4x4 jacobian, denoted j4
    keep_lattice = kwargs.pop('keep_lattice', False)
    convergence = kwargs.pop('convergence', DConstant.OrbConvergence)
    max_iterations = kwargs.pop('max_iterations', DConstant.OrbMaxIter)
    xy_step = kwargs.pop('XYStep', DConstant.XYStep)
    ref_in = numpy.zeros((6,)) if guess is None else numpy.copy(guess)
    ref_in[4] = 0.0 if dp is None else dp

    scaling = xy_step * numpy.array([1.0, 1.0, 1.0, 1.0])
    delta_matrix = numpy.zeros((6, 5), order='F')
    for i in range(4):
        delta_matrix[i, i] = scaling[i]
    id4 = numpy.asfortranarray(numpy.identity(4))
    change = 1
    itercount = 0
    while (change > convergence) and itercount < max_iterations:
        in_mat = ref_in.reshape((6, 1)) + delta_matrix
        _ = lattice_pass(ring, in_mat, refpts=[], keep_lattice=keep_lattice)
        # the reference particle after one turn
        ref_out = in_mat[:, 4]
        # 4x4 jacobian matrix from numerical differentiation:
        # f(x+d) - f(x) / d
        j4 = (in_mat[:4, :4] - in_mat[:4, 4:]) / scaling
        a = j4 - id4  # f'(r_n) - 1
        b = ref_out[:4] - ref_in[:4]
        b_over_a = numpy.linalg.solve(a, b)
        r_next = ref_in - numpy.append(b_over_a, numpy.zeros((2,)))
        # determine if we are close enough
        change = numpy.linalg.norm(r_next - ref_in)
        itercount += 1
        ref_in = r_next
        keep_lattice = True

    if itercount == max_iterations:
        warnings.warn(AtWarning('Maximum number of iterations reached. '
                                'Possible non-convergence'))
    return ref_in


@check_radiation(False)
def _orbit_dct(ring, dct=None, guess=None, **kwargs):
    """Solver for fixed path lengthening"""
    keep_lattice = kwargs.pop('keep_lattice', False)
    convergence = kwargs.pop('convergence', DConstant.OrbConvergence)
    max_iterations = kwargs.pop('max_iterations', DConstant.OrbMaxIter)
    xy_step = kwargs.pop('XYStep', DConstant.XYStep)
    ref_in = numpy.zeros((6,)) if guess is None else numpy.copy(guess)

    scaling = xy_step * numpy.array([1.0, 1.0, 1.0, 1.0, 1.0])
    delta_matrix = numpy.zeros((6, 6), order='F')
    for i in range(5):
        delta_matrix[i, i] = scaling[i]
    theta5 = numpy.zeros((5,), order='F')
    theta5[4] = 0.0 if dct is None else dct
    id5 = numpy.zeros((5, 5), order='F')
    for i in range(4):
        id5[i, i] = 1.0
    idx = numpy.array([0, 1, 2, 3, 5])
    change = 1
    itercount = 0
    while (change > convergence) and itercount < max_iterations:
        in_mat = ref_in.reshape((6, 1)) + delta_matrix
        _ = lattice_pass(ring, in_mat, refpts=[], keep_lattice=keep_lattice)
        # the reference particle after one turn
        ref_out = in_mat[:, -1]
        # 5x5 jacobian matrix from numerical differentiation:
        # f(x+d) - f(x) / d
        j5 = (in_mat[idx, :5] - in_mat[idx, 5:]) / scaling
        a = j5 - id5  # f'(r_n) - 1
        b = ref_out[idx] - numpy.append(ref_in[:4], 0.0) - theta5
        b_over_a = numpy.linalg.solve(a, b)
        r_next = ref_in - numpy.append(b_over_a, 0.0)
        # determine if we are close enough
        change = numpy.linalg.norm(r_next - ref_in)
        itercount += 1
        ref_in = r_next
        keep_lattice = True

    if itercount == max_iterations:
        warnings.warn(AtWarning('Maximum number of iterations reached. '
                                'Possible non-convergence'))
    return ref_in


def find_orbit4(ring, dp=0.0, refpts=None, dct=None, orbit=None,
                keep_lattice=False, **kwargs):
    """findorbit4 finds the closed orbit in the 4-d transverse phase
    space by numerically solving for a fixed point of the one turn
    map M calculated with lattice_pass.

        (X, PX, Y, PY, dP, CT2 ) = M (X, PX, Y, PY, dP, CT1)

    under the CONSTANT MOMENTUM constraint dP and with NO constraint
    on the 6-th coordinate CT

    IMPORTANT!!! findorbit4 imposes a constraint on dP and relaxes
    the constraint on the revolution frequency. A physical storage
    ring does exactly the opposite: the momentum deviation of a
    particle on the closed orbit settles at the value
    such that the revolution is synchronous with the RF cavity

                HarmNumber*Frev = Frf

    To impose this artificial constraint in find_orbit4, PassMethod
    used for any element SHOULD NOT
    1.  change the longitudinal momentum dP (cavities , magnets with radiation)
    2.  have any time dependence (localized impedance, fast kickers etc)

    PARAMETERS
        ring            lattice description (radiation must be OFF)
        dp              momentum deviation. Defaults to 0
        refpts          elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.
                        Defaults to None, if refpts is None an empty array is
                        returned for orbit.

    OUTPUT
        orbit0          ((6,) closed orbit vector at the entrance of the
                        1-st element (x,px,y,py)
        orbit           (6, Nrefs) closed orbit vector at each location
                        specified in refpts

    KEYWORDS
        dct=None        path lengthening. If specified, dp is ignored and
                        the off-momentum is deduced from the path lengthening.
        orbit=None      avoids looking for initial the closed orbit if is
                        already known ((6,) array). find_orbit4 propagates it
                        to the specified refpts.
        guess           (6,) initial value for the closed orbit. It may help
                        convergence. Default: (0, 0, 0, 0, 0, 0)
        keep_lattice    Assume no lattice change since the previous tracking.
                        Default: False
        convergence     Convergence criterion. Default: 1.e-12
        max_iterations  Maximum number of iterations. Default: 20
        XYStep          Step size. Default: DConstant.XYStep

    See also find_sync_orbit, find_orbit6.
    """
    if orbit is None:
        if dct is not None:
            orbit = _orbit_dct(ring, dct, keep_lattice=keep_lattice, **kwargs)
        else:
            orbit = _orbit_dp(ring, dp, keep_lattice=keep_lattice, **kwargs)
        keep_lattice = True

    uint32refs = uint32_refpts(refpts, len(ring))
    # bug in numpy < 1.13
    all_points = numpy.empty((0, 6), dtype=float) if len(
        uint32refs) == 0 else numpy.squeeze(
        lattice_pass(ring, orbit.copy(order='K'), refpts=uint32refs,
                     keep_lattice=keep_lattice), axis=(1, 3)).T
    return orbit, all_points


def find_sync_orbit(ring, dct=0.0, refpts=None, dp=None, orbit=None,
                    keep_lattice=False, **kwargs):
    """find_sync_orbit finds the closed orbit, synchronous with the RF cavity
    and momentum deviation dP (first 5 components of the phase space vector)
    % by numerically solving  for a fixed point
    % of the one turn map M calculated with lattice_pass

        (X, PX, Y, PY, dP, CT2 ) = M (X, PX, Y, PY, dP, CT1)

    under the constraint dCT = CT2 - CT1 = C/Frev - C/Frev0, where
    Frev0 = Frf0/HarmNumber is the design revolution frequency
    Frev  = (Frf0 + dFrf)/HarmNumber is the imposed revolution frequency

    IMPORTANT!!!  find_sync_orbit imposes a constraint (CT2 - CT1) and
    dP2 = dP1 but no constraint on the value of dP1, dP2
    The algorithm assumes time-independent fixed-momentum ring
    to reduce the dimensionality of the problem.

    To impose this artificial constraint in find_sync_orbit
    PassMethod used for any element SHOULD NOT
    1.  change the longitudinal momentum dP (cavities , magnets with radiation)
    2.  have any time dependence (localized impedance, fast kickers etc).

    PARAMETERS
        ring            lattice description (radiation must be OFF)
        dct             Path length deviation. Default: 0
        refpts          elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.
                        Defaults to None, if refpts is None an empty array is
                        returned for orbit.

    OUTPUT
        orbit0          ((6,) closed orbit vector at the entrance of the
                        1-st element (x,px,y,py)
        orbit           (6, Nrefs) closed orbit vector at each location
                        specified in refpts

    KEYWORDS
        orbit=None      avoids looking for initial the closed orbit if is
                        already known ((6,) array). find_sync_orbit propagates
                        it to the specified refpts.
        guess           (6,) initial value for the closed orbit. It may help
                        convergence. Default: (0, 0, 0, 0, 0, 0)
        keep_lattice    Assume no lattice change since the previous tracking.
                        Default: False
        convergence     Convergence criterion. Default: 1.e-12
        max_iterations  Maximum number of iterations. Default: 20
        XYStep          Step size. Default: DConstant.XYStep

    See also find_orbit4, find_orbit6.
    """
    if orbit is None:
        if dp is not None:
            orbit = _orbit_dp(ring, dp, keep_lattice=keep_lattice, **kwargs)
        else:
            orbit = _orbit_dct(ring, dct,  keep_lattice=keep_lattice, **kwargs)
        keep_lattice = True

    uint32refs = uint32_refpts(refpts, len(ring))
    # bug in numpy < 1.13
    all_points = numpy.empty((0, 6), dtype=float) if len(
        uint32refs) == 0 else numpy.squeeze(
        lattice_pass(ring, orbit.copy(order='K'), refpts=uint32refs,
                     keep_lattice=keep_lattice), axis=(1, 3)).T
    return orbit, all_points


def _orbit6(ring, cavpts=None, guess=None, keep_lattice=False, **kwargs):
    """Solver for 6D motion"""
    convergence = kwargs.pop('convergence', DConstant.OrbConvergence)
    max_iterations = kwargs.pop('max_iterations', DConstant.OrbMaxIter)
    xy_step = kwargs.pop('XYStep', DConstant.XYStep)
    dp_step = kwargs.pop('DPStep', DConstant.DPStep)
    method = kwargs.pop('method', ELossMethod.TRACKING)

    l0 = get_s_pos(ring, len(ring))[0]
    f_rf = ring.get_rf_frequency()
    harm_number = round(f_rf*l0/clight)

    if guess is None:
        _, dt = get_timelag_fromU0(ring, method=method, cavpts=cavpts)
        # Getting timelag by tracking uses a different lattice,
        # so we cannot now use the same one again.
        if method is ELossMethod.TRACKING:
            keep_lattice = False
        ref_in = numpy.zeros((6,), order='F')
        ref_in[5] = -dt
    else:
        ref_in = numpy.copy(guess)

    theta = numpy.zeros((6,))
    theta[5] = clight * harm_number / f_rf - l0

    scaling = xy_step * numpy.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0]) + \
              dp_step * numpy.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
    delta_matrix = numpy.asfortranarray(
        numpy.concatenate((numpy.diag(scaling), numpy.zeros((6, 1))), axis=1))

    id6 = numpy.asfortranarray(numpy.identity(6))
    change = 1
    itercount = 0
    while (change > convergence) and itercount < max_iterations:
        in_mat = ref_in.reshape((6, 1)) + delta_matrix
        _ = lattice_pass(ring, in_mat, refpts=[], keep_lattice=keep_lattice)
        # the reference particle after one turn
        ref_out = in_mat[:, 6]
        # 6x6 jacobian matrix from numerical differentiation:
        # f(x+d) - f(x) / d
        j6 = (in_mat[:, :6] - in_mat[:, 6:]) / scaling
        a = j6 - id6  # f'(r_n) - 1
        b = ref_out[:] - ref_in[:] - theta
        # b_over_a, _, _, _ = numpy.linalg.lstsq(a, b, rcond=-1)
        b_over_a = numpy.linalg.solve(a, b)
        r_next = ref_in - b_over_a
        # determine if we are close enough
        change = numpy.linalg.norm(r_next - ref_in)
        itercount += 1
        ref_in = r_next
        keep_lattice = True

    if itercount == max_iterations:
        warnings.warn(AtWarning('Maximum number of iterations reached. '
                                'Possible non-convergence'))
    return ref_in


def find_orbit6(ring, refpts=None, orbit=None, dp=None, dct=None,
                keep_lattice=False, **kwargs):
    """find_orbit6 finds the closed orbit in the full 6-D phase space
    by numerically solving  for a fixed point of the one turn
    map M calculated with lattice_pass

    (X, PX, Y, PY, DP, CT2 ) = M (X, PX, Y, PY, DP, CT1)

    with constraint  CT2 - CT1 = C*HarmNumber(1/Frf - 1/Frf0)

    IMPORTANT!!! find_orbit6 is a realistic simulation
    1.  The Frf frequency in the RF cavities (may be different from Frf0)
        imposes the synchronous condition
        CT2 - CT1 = C*HarmNumber(1/Frf - 1/Frf0)
    2.  The algorithm numerically calculates
        6-by-6 Jacobian matrix J6. In order for (J-E) matrix
        to be non-singular it is NECESSARY to use a realistic
        PassMethod for cavities with non-zero momentum kick
        (such as ThinCavityPass).
    3.  find_orbit6 can find orbits with radiation.
        In order for the solution to exist the cavity must supply
        adequate energy compensation.
        In the simplest case of a single cavity, it must have
        'Voltage' field set so that Voltage > Erad - energy loss per turn
    4.  There is a family of solutions that correspond to different RF buckets
        They differ in the 6-th coordinate by C*Nb/Frf. Nb = 1 .. HarmNum-1
    5.  The value of the 6-th coordinate found at the cavity gives
        the equilibrium RF phase. If there is no radiation the phase is 0;

    PARAMETERS
        ring            lattice description (radiation must be ON)
        refpts          elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.
                        Defaults to None, if refpts is None an empty array is
                        returned for orbit.

    OUTPUT
        orbit0          ((6,) closed orbit vector at the entrance of the
                        1-st element (x,px,y,py)
        orbit           (6, Nrefs) closed orbit vector at each location
                        specified in refpts

    KEYWORDS
        orbit=None      avoids looking for initial the closed orbit if is
                        already known ((6,) array). find_orbit6 propagates it
                        to the specified refpts.
        guess           Initial value for the closed orbit. It may help
                        convergence. The default is computed from the energy
                        loss of the ring
        keep_lattice    Assume no lattice change since the previous tracking.
                        Default: False
        method          Method for energy loss computation (see get_energy_loss)
                        default: ELossMethod.TRACKING
        cavpts=None     Cavity location. If None, use all cavities. This is used
                        to compute the initial synchronous phase.
        convergence     Convergence criterion. Default: 1.e-12
        max_iterations  Maximum number of iterations. Default: 20
        XYStep          Step size. Default: DConstant.XYStep
        DPStep          Step size. Default: DConstant.DPStep

    See also find_orbit4, find_sync_orbit.
    """
    if not (dp is None and dct is None):
        warnings.warn(AtWarning('In 6D, "dp" and "dct" are ignored'))
    if orbit is None:
        orbit = _orbit6(ring, keep_lattice=keep_lattice, **kwargs)
        keep_lattice = True

    uint32refs = uint32_refpts(refpts, len(ring))
    # bug in numpy < 1.13
    all_points = numpy.empty((0, 6), dtype=float) if len(
        uint32refs) == 0 else numpy.squeeze(
        lattice_pass(ring, orbit.copy(order='K'), refpts=uint32refs,
                     keep_lattice=keep_lattice), axis=(1, 3)).T
    return orbit, all_points


def find_orbit(ring, refpts=None, **kwargs):
    """find_orbit finds the closed orbit by numerically getting the fixed point
    of the one turn map M calculated with lattice_pass.

    Depending on the the lattice, find_orbit will:
    - use find_orbit6 if ring.radiation is ON,
    - use find_sync_orbit if ring.radiation is OFF and dct is specified,
    - use find_orbit4 otherwise

    PARAMETERS
        ring            Sequence of AT elements
        refpts          elements at which data is returned.

    OUTPUT
        orbit0          ((6,) closed orbit vector at the entrance of the
                        1-st element
        orbit           (6, Nrefs) closed orbit vector at each location
                        specified in refpts

    KEYWORDS
        dp=0            Momentum deviation, when radiation is OFF
        dct=0            Path lengthening, when radiation ids OFF
        keep_lattice    Assume no lattice change since the previous tracking.
                        Default: False
        guess=None      Initial guess for the closed orbit. It may help
                        convergence.
        orbit=None      Orbit at entrance of the lattice, if known. find_orbit
                        will then propagate it to the selected reference points
        For other keywords, refer to the underlying methods

    See also find_orbit4, find_sync_orbit, find_orbit6
    """
    if ring.radiation:
        return find_orbit6(ring, refpts=refpts, **kwargs)
    else:
        return find_orbit4(ring, refpts=refpts, **kwargs)


Lattice.find_orbit4 = find_orbit4
Lattice.find_sync_orbit = find_sync_orbit
Lattice.find_orbit6 = find_orbit6
Lattice.find_orbit = find_orbit
