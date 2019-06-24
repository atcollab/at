import os
import numpy
import matlab
import at
from at.physics.rad_int import get_synchrotron_integrals
from at.physics.rad_int import calc_disp
from at.physics.rad_int import calc_twiss
from at.physics.rad_int import calc_twiss1
from at.physics.rad_int import calc_mx
from at.physics.rad_int import calc_rad_int


def _test_find_radiation_integrals(engine):
    path = os.path.join(os.path.dirname(__file__),
                        '../../../atip/rings/diad.mat')
    py_lattice = at.load_mat(path)
    ml_lattice = engine.load(path)['RING']
    twiss_data = at.physics.get_twiss(py_lattice, refpts=range(len(py_lattice)),
                                      get_chrom=True, ddp=1e-5)
    """"""
    for i in range(len(py_lattice)):
        elem = py_lattice[i]
        if isinstance(elem, at.lattice.elements.Dipole):
            theta = elem.BendingAngle
            rho = elem.Length / theta
            K = elem.K
            a0 = twiss_data[3]['alpha'][i, 0]
            b0 = twiss_data[3]['beta'][i, 0]
            D0 = twiss_data[3]['dispersion'][i, 0]
            D0p = twiss_data[3]['dispersion'][i, 1] + (numpy.tan(elem.EntranceAngle)*D0)/rho
            #py_data = calc_disp(rho, theta, D0, D0p, K)
            #ml_data = numpy.squeeze(engine.calcdisp(rho, theta, D0, D0p, K))
            #py_data = calc_twiss(rho, theta, a0, b0, K)
            #ml_data = numpy.squeeze(engine.calctwiss(rho, theta, a0, b0, K))
            py_data = calc_rad_int(elem, a0, b0, twiss_data[3]['dispersion'][i])
            ml_data = numpy.squeeze(engine.calcRadInt(rho, theta, a0, b0, D0, D0p,K, elem.EntranceAngle, elem.ExitAngle))
            try:
                numpy.testing.assert_almost_equal(py_data, ml_data, decimal=10)
            except AssertionError as e:
                raise e
                #raise KeyError(i, e)
    raise KeyError('Sucessful')#\npydata{}\nmldata{}'.format(py_data, ml_data))
    """"""
    """Full integrals"""
    # Python call
    py_integrals = get_synchrotron_integrals(py_lattice)
    # Matlab call
    results = engine.atsummary(ml_lattice)
    ml_integrals = numpy.squeeze(results['integrals'])
    # Comparison
    numpy.testing.assert_almost_equal(py_integrals, ml_integrals[:5])
