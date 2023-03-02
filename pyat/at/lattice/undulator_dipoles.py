import at
import copy


def undulator_dipoles(L, nperiod,
              BendAngle=0.0,
              B0=None,
              Energy=None,
              magnetmodel='rectangularbend',
              PoleDistance=0.0, name=''):
    """

    :param L:
    :param nperiod:
    :param BendAngle:
    :param B0:
    :param Energy:
    :param magnetmodel:
    :param PoleDistance:
    :return:

    input:
       Lund= undulator length
       nperiod = number of periods
       'BendAngle', value : half pole bending angle in rad
       'B0andEnergy', value (2x1): [half pole B0 field in T, Energy in eV]
                                   converts to bending angle in rad.
       'magnetmodel', value : 'multipoles' (default) or 'rectangularbend'
       'PoleGap', value : drift space between gaps defualt (0.0)

     if neither BendAngle nor B0andEnergy are provided then 'BendAngle' is 0.0

     output:
        cellarray of elements describing un undulator of length LUnd divided
        in nperiod periods, each described as follows:
            [negpole,drift,pospole,pospole,drift,negpole]
        if 'PoleGap' is 0.0 (default), then
        [negpole,pospole,pospole,negpole]

    example:
    1)  und=atundulator(1.6,61,'B0andEnergy',[0.4 6.04e9])
    2)  und=atundulator(1.6,61,'BendAngle',-0.007984472464733)
    3)  und=atundulator(1.6,61,'B0andEnergy',[0.4 6.04e9],'magnetmodel','rectangularbend')
    4)  und=atundulator(1.6,61,'B0andEnergy',[0.4 6.04e9],'PoleGap',0.001);
    """

    # length of one period
    periodL = L / nperiod

    DistPole = PoleDistance
    LPole = (periodL - 2 * DistPole) / 4

    if (B0 is not None) and (Energy is not None):
        B = B0
        Brho = Energy / 299792458
        AngPole = (B * LPole) / Brho
    elif BendAngle is not None:
        AngPole = BendAngle
    else:
        AngPole = 0.0
    print(f'half PoleAngle = {AngPole}')

    # function to define one period of the undulator
    def makeundperiod(halfnegpole, halfpospole, driftpole):

        if driftpole.Length > 0:
            undper = [
                halfnegpole.deepcopy(),
                driftpole.deepcopy(),
                halfpospole.deepcopy(),
                halfpospole.deepcopy(),
                driftpole.deepcopy(),
                halfnegpole.deepcopy(),
                ]
        elif driftpole.Length == 0:
            undper = [
                halfnegpole.deepcopy(),
                halfpospole.deepcopy(),
                halfpospole.deepcopy(),
                halfnegpole.deepcopy(),
                ]

        return undper

    if magnetmodel == 'multipoles':

        undperiod = makeundperiod(
            at.Multipole('NegPole', LPole, -AngPole / LPole),
            at.Multipole('PosPole', LPole,  AngPole / LPole),
            at.Drift('PoleGap', DistPole))

    elif magnetmodel == 'rectangularbend':

        undperiod = makeundperiod(
            at.Dipole('NegPole', LPole, bending_angle=-AngPole,
                      EntranceAngle=-AngPole/2, ExitAngle=-AngPole/2),
            at.Dipole('PosPole', LPole, bending_angle=AngPole,
                      EntranceAngle=AngPole/2, ExitAngle=AngPole/2),
            at.Drift('PoleGap', DistPole))

    # repeat undperiod n period times and amke all elements independent
    undulator = at.Lattice([el.deepcopy() for el in undperiod * nperiod], periodicity=1, energy=Energy)

    # place a marker at the center of the undulator
    undulator.sbreak(L/2, break_elems=at.Marker(f'{name}_center'))

    return undulator


def insert_undulator(ring, und, index, matchoptics=False):
    """
    insert an undulator created using undulator.py at the specified location
    
    :param ring: AT lattice where to place the ring
    :param und: AT lattice structure obtained by undulator_dipoles
    :param index: location at wich the undulator is placed. Half length will be taken upstream and half downstream
    :param matchoptics: use quadrupoles in positive and negavite poles of undulator to match optics
    :return: ring with undulator and matched quadrupoles, ring with undulators 
    """

    if type(ring[index]) is not at.elements.Marker:
        raise Exception('location to place undulator is not a marker')
    
    und = at.Lattice(und, energy=ring.energy, periodicity=1)

    ring_und = ring.deepcopy()
    # look if there is space
    L_before = ring_und[index - 1].Length
    L_after = ring_und[index + 1].Length
    if L_before > und.circumference / 2 and  L_after > und.circumference / 2:
        # insert undulator
        # substract half length on each side
        ring_und[index - 1].Length = L_before - und.circumference / 2
        ring_und[index + 1].Length = L_after - und.circumference / 2
        ring_und = at.Lattice(ring_und[0:index] + und.deepcopy() + ring_und[index:], periodicity=1)

    else:
        raise Exception(f'not enough space at this location: {L_before} + {L_after} meters')

    ring_und_matched = ring_und

    if matchoptics:
        # match quadrupoles
        # Define the variables, one quad at entrance and one at exit of undulator
        indp = at.get_refpts(ring_und, 'PosPole*')[0]  #[[0, -1]]
        indn = at.get_refpts(ring_und, 'NegPole*')[-1]  #[[0, -1]]
    
        # add some PolynomB in case changes relative to 0.
        for ind, el in enumerate(ring_und):
            if el.FamName == 'PosPole':
                el.MaxOrder = 1
                el.PolynomB[1] = 0.0
            if el.FamName == 'NegPole':
                el.MaxOrder = 1
                el.PolynomB[1] = 0.0
        
        variables = [
            at.ElementVariable(indp, 'PolynomB', index=1, name='focussing'),
            at.ElementVariable(indn, 'PolynomB', index=1, name='defocussing')]
    
        match_pos = len(ring_und)  #
        twiss_in, _, _ = at.linopt4(ring, 0)
    
        cst1 = at.LinoptConstraints(ring_und, twiss_in=twiss_in)
    
        l0, _, _ = at.linopt4(ring, len(ring))
        cst1.add('alpha', l0.alpha, refpts=match_pos, name='alpha_end', weight=1e-3)
    
        ring_und_matched = at.match(ring_und_matched, variables, (cst1,), method='trf', diff_step=1e-10)  # diff_step=1e-5
        # ring_und_matched = at.match(ring_und_matched, variables, (cst1,), method='lm')  # diff_step

    return ring_und_matched, ring_und