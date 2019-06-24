"""
Synchrotron integrals
Taken as a mixure of atsummary.m and ringpara.m
"""
import numpy
from numpy import tan, cos, sin, cosh, sinh
from at.physics import get_twiss
from at.lattice import elements


def get_synchrotron_integrals(ring, twiss=None):
    I1 = 0.0
    I2 = 0.0
    I3 = 0.0
    I4 = 0.0
    I5 = 0.0
    if twiss is None:
        twiss = get_twiss(ring, refpts=range(len(ring)), get_chrom=True,
                          ddp=1e-5)
    for i in range(len(ring)):
        elem = ring[i]
        if isinstance(elem, elements.Dipole):
            if elem.BendingAngle != 0:
                int_contrib = calc_rad_int(elem, twiss[3]['alpha'][i, 0],
                                           twiss[3]['beta'][i, 0],
                                           twiss[3]['dispersion'][i])
                I1 += int_contrib[0]  # dI1
                I2 += int_contrib[1]  # dI2
                I3 += int_contrib[2]  # dI3
                I4 += int_contrib[3]  # dI4
                I5 += int_contrib[4]  # dI5
    return I1, I2, I3, I4, I5


def calc_rad_int(elem, alpha, beta, dispersion):
    # Calcuate a dipole's contribution to the radiation integrals
    rho = elem.Length / elem.BendingAngle  # radius of the dipole
    th1 = elem.EntranceAngle
    th2 = elem.ExitAngle
    # Entrance edge focusing
    D0p = dispersion[1] + (tan(th1) * dispersion[0])/rho
    a0 = alpha - (tan(th1) * beta)/rho
    # Split the dipole into N pieces
    n = 100
    th = numpy.array(range(n+1)) * (elem.BendingAngle/n)
    Dx = numpy.zeros(n+1)
    Dxp = numpy.zeros(n+1)
    curH = numpy.zeros(n+1)
    # Compute Twiss parameters inside the dipole
    for i in range(n+1):
        Dx[i], Dxp[i] = calc_disp(rho, th[i], dispersion[0], D0p, elem.K)
        ax, bx = calc_twiss(rho, th[i], a0, beta, elem.K)
        curH[i] = (Dx[i]**2 + (ax*Dx[i] + bx*Dxp[i])**2)/bx
    # Exit edge focusing
    Dxp[-1] = Dxp[-1] + (tan(th2)*Dx[-1])/rho
    ax = ax - (tan(th2)*bx)/rho
    curH[-1] = (Dx[-1]**2 + (ax*Dx[-1] + bx*Dxp[-1])**2)/bx
    # Average data
    curHavg = ((curH[0] + curH[-1])/2 + sum(curH[1:-2]))/n
    Dxavg = ((Dx[0] + Dx[-1])/2 + sum(Dx[1:-2]))/n
    # Integral contributions
    dI1 = Dxavg * elem.BendingAngle
    dI2 = abs(elem.BendingAngle / rho)
    dI3 = abs(elem.BendingAngle / rho**2)
    dI4 = dI1*(2*elem.K + rho**-2) - (tan(th1)*Dx[0])/rho**2 - (tan(th2)*Dx[-1])/rho**2
    dI5 = curHavg * dI3
    return Dx
    # return curH[0], curH[-1], sum(curH[1:-2]), n
    return dI1, dI2, dI3, dI4, dI5, curHavg, Dxavg


def calc_disp(rho, theta, D0, D0p, K):
    # calcualte dispersion function inside a combined-function dipole
    s = rho * theta
    if K > (-rho**-2):  # horizontal focusing
        sqK = numpy.sqrt(K + rho**-2)
        Dx = D0*cos(sqK*s) + (sin(sqK*s)*D0p)/sqK + ((1 - cos(sqK*s))/rho)/sqK**2
        Dxp = -D0*sqK*sin(sqK*s) + D0p*cos(sqK*s) + (sin(sqK*s)/rho)/sqK
    else:  # horizontal defocusing
        sqK = numpy.sqrt(-(K + rho**-2))
        Dx = D0*cosh(sqK*s) + (sinh(sqK*s)*D0p)/sqK + ((cosh(sqK*s) - 1)/rho)/sqK**2
        Dxp = D0*sqK*sinh(sqK*s) + D0p*cosh(sqK*s) + (sinh(sqK*s)/rho)/sqK
    return Dx, Dxp


def calc_twiss(rho, theta, a0, b0, K):
    # twx2 not used as calculation of ax and bx done directly from Mx, a, b, g
    # calculate twiss function inside a combined-function dipole manget
    s = rho * theta
    g0 = (1 + a0**2)/b0
    if K > (-rho**-2):  # horizontal focusing
        sqK = numpy.sqrt(K + rho**-2)
        ax = b0*cos(sqK*s)*sqK*sin(sqK*s) + a0*cos(sqK*s)**2 - a0*sin(sqK*s)**2 - (g0*sin(sqK*s)*cos(sqK*s))/sqK
        bx = b0*cos(sqK*s)**2 - (a0*2*cos(sqK*s)*sin(sqK*s))/sqK + g0*(sin(sqK*s)/sqK)**2
    else:  # horizontal defocusing
        sqK = numpy.sqrt(-(K + rho**-2))
        ax = -b0*cosh(sqK*s)*sqK*sinh(sqK*s) + a0*cosh(sqK*s)**2 + a0*sinh(sqK*s)**2 - (g0*sinh(sqK*s)*cosh(sqK*s))/sqK
        bx = b0*cosh(sqK*s)**2 - (a0*2*cosh(sqK*s)*sinh(sqK*s))/sqK + g0*(sinh(sqK*s)/sqK)**2
    return ax, bx


def calc_twiss1(rho, theta, a0, b0, K):
    # calculate twiss function inside a combined-function dipole manget
    s = rho * theta
    if K > (-rho**-2):  # horizontal focusing
        sqK = numpy.sqrt(K + rho**-2)
        Mx = [cos(sqK*s), sin(sqK*s)/sqK, -sqK*sin(sqK*s), cos(sqK*s)]
    else:  # horizontal defocusing
        sqK = numpy.sqrt(-(K + rho**-2))
        Mx = [cosh(sqK*s), sinh(sqK*s)/sqK, sqK*sinh(sqK*s), cosh(sqK*s)]
    g0 = (1 + a0**2)/b0
    twx2 = numpy.array([[Mx[0]**2, -2*Mx[0]*Mx[1], Mx[1]**2],
                        [-Mx[0]*Mx[2], Mx[0]*Mx[3] + Mx[1]*Mx[2], -Mx[1]*Mx[3]],
                        [Mx[2]**2, -2*Mx[2]*Mx[3], Mx[3]**2]]).dot([b0, a0, g0])
    # ax = -Mx[0]*Mx[2]*b0 + Mx[0]*Mx[3]*a0 + Mx[1]*Mx[2]*a0 - Mx[1]*Mx[3]*g0
    # bx = b0*Mx[0]**2 - 2*Mx[0]*Mx[1]*a0 + g0*Mx[1]**2
    return twx2[1], twx2[0]


def calc_mx(rho, K, theta):
    # redundant as directly included in both calc twiss functions
    s = rho * theta
    if K > (-rho**-2):  # horizontal focusing
        sqK = numpy.sqrt(K + rho**-2)
        Mx = [cos(sqK*s), sin(sqK*s)/sqK, -sqK*sin(sqK*s), cos(sqK*s)]
    else:  # horizontal defocusing
        sqK = numpy.sqrt(-(K + rho**-2))
        Mx = [cosh(sqK*s), sinh(sqK*s)/sqK, sqK*sinh(sqK*s), cosh(sqK*s)]
    return Mx
