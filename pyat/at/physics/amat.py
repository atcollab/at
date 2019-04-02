""""""
from collections import namedtuple
import numpy
from scipy.linalg import block_diag, eig, inv, det
from math import pi
from ..physics import jmat

_i2 = numpy.array([[-1.j, -1.], [1., 1.j]])
_vxyz = [_i2, block_diag(_i2, _i2), block_diag(_i2, _i2, _i2)]
_submat = [slice(0, 2), slice(2, 4), slice(6, 3, -1)]

_Data1 = namedtuple('R66Data', ('tunes', 'damping_rates', 'mode_matrices'))
_Data2 = namedtuple('R66Data', ('tunes', 'damping_rates', 'mode_matrices',
                                'mode_emittances'))


def amat(tt):
    """
    find A matrix from one turn map matrix T such that:

               [Rotx  0    0  ]
    inv(A)*T*A=[ 0   Rotz  0  ]
               [ 0    0   Rots]

    Order it so that it is close to the order of x,y,z
    also ensure that positive modes are before negative so that
    one has proper symplecticity
    B. Nash July 18, 2013

    INPUT
    T       (m, m)  transfer matrix for 1 turn

    OUTPUT
    A       (m, m)  A-matrix

    """
    nv = tt.shape[0]
    dms = int(nv / 2)
    jmt = jmat(dms)
    select = numpy.arange(0, nv, 2)
    rbase = numpy.stack((select, select), axis=1).flatten()
    # noinspection PyTupleAssignmentBalance
    _, vv = eig(tt)
    # Compute the norms
    vp = numpy.dot(vv.conj().T, jmt)
    n = -0.5j * numpy.sum(vp.T * vv, axis=0)
    # Move positive before negatives
    order = rbase + (n < 0)
    vv = vv[:, order]
    n = n[order]
    # Normalize vectors
    vn = vv / numpy.sqrt(abs(n)).reshape((1, nv))
    # find the vectors that project most onto x,y,z, and reorder
    # nn will have structure
    #  n1x n1y n1z
    #  n2x n2y n2z
    #  n3x n3y n3z
    nn = 0.5 * abs(numpy.sqrt(-1.j * vn.conj().T.dot(jmt).dot(_vxyz[dms - 1])))
    ind = numpy.argmax(nn[select, :][:, select], axis=0)
    v_ordered = vn[:, 2 * ind]
    aa = numpy.vstack((numpy.real(v_ordered), numpy.imag(v_ordered))).reshape(
        (nv, nv), order='F')
    return aa


def get_mode_matrices(a):
    """Given a (m, m) A matrix , find the m normal modes"""
    dms = int(a.shape[0] / 2)
    return numpy.stack([numpy.dot(a[:, s], a.T[s, :]) for s in _submat[:dms]],
                       axis=0)


def get_tunes_damp(tt, rr):
    """
    mode_emit, damping_rates, tunes = get_tunes_damp(T, R)

    INPUT
        T                   (m, m) transfer matrix for 1 turn
        R                   (m, m) beam matrix (optional)

        m can be 2 (single plane), 4 (betatron motion) or 6 (full motion)

    OUTPUT
        named tuple with the follwing attributes:
        tunes               (m,) tunes of the m normal modes
        damping_rates       (m,) damping rates of the m normal modes
        mode_matrices       (3, m, m) the R-matrices of the m normal modes
        mode_emittances     Only if R is specified: (m,) emittance of each mode
    """

    def decode(rot22):
        tune = (numpy.arctan2(rot22[0, 1] - rot22[1, 0],
                              rot22[0, 0] + rot22[1, 1]) / 2.0 / pi) % 1
        chi = -numpy.log(numpy.sqrt(det(rot22)))
        return chi, tune

    nv = tt.shape[0]
    dms = int(nv / 2)
    jmt = jmat(dms)
    aa = amat(tt)
    rmat = inv(aa).dot(tt.dot(aa))
    damping_rates, tunes = zip(*(decode(rmat[s, s]) for s in _submat[:dms]))
    if rr is None:
        return _Data1(tunes, damping_rates, get_mode_matrices(aa))
    else:
        rdiag = numpy.diag(aa.T.dot(jmt.dot(rr.dot(jmt.dot(aa)))))
        mode_emit = -0.5 * (rdiag[0:nv:2] + rdiag[1:nv:2])
        return _Data2(tunes, damping_rates, get_mode_matrices(aa), mode_emit)
