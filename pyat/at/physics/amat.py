""""""
import numpy
from scipy.linalg import block_diag, eig, inv, solve
from math import pi
from at.lattice import AtError

__all__ = ['a_matrix', 'amat', 'jmat', 'jmatswap',
           'get_tunes_damp', 'get_mode_matrices', 'symplectify']

_i2 = numpy.array([[-1.j, -1.], [1., 1.j]])

# Prepare symplectic identity matrix
_j2 = numpy.array([[0., 1.], [-1., 0.]])
_jm = [_j2, block_diag(_j2, _j2), block_diag(_j2, _j2, _j2)]
_jmswap = [_j2, block_diag(_j2, _j2), block_diag(_j2, _j2, _j2.T)]

_vxyz = [_i2, block_diag(_i2, _i2), block_diag(_i2, _i2, _i2)]
_submat = [slice(0, 2), slice(2, 4), slice(4, 6)]


def jmat(ind):
    """
    Return the antisymetric block diagonal matrix [[0, 1][-1, 0]]

    INPUT
        ind     1, 2 or 3. Matrix dimension

    OUTPUT
        jm      block diagonal matrix, (2, 2) or (4, 4) or (6, 6)
    """
    return _jm[ind - 1]


def jmatswap(ind):
    """Modified version of jmat to deal with the swap of the
    longitudinal coordinates"""
    return _jmswap[ind - 1]


def a_matrix(tt):
    """A, eigval = a_matrix(T)
    Find the A matrix from one turn map matrix T such that:

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
    eigval  (m,)    vector of Eigen values of T
    """
    nv = tt.shape[0]
    dms = int(nv / 2)
    jmt = jmatswap(dms)
    select = numpy.arange(0, nv, 2)
    rbase = numpy.stack((select, select), axis=1).flatten()

    # noinspection PyTupleAssignmentBalance
    lmbd, vv = eig(tt)
    # Compute the norms
    vp = numpy.dot(vv.conj().T, jmt)
    n = -0.5j * numpy.sum(vp.T * vv, axis=0)
    if any(abs(n) < 1.0E-12):
        raise AtError('Unstable ring')
    # Move positive before negatives
    order = rbase + (n < 0)
    vv = vv[:, order]
    n = n[order]
    lmbd = lmbd[order]
    # Normalize vectors
    vn = vv / numpy.sqrt(abs(n)).reshape((1, nv))
    # find the vectors that project most onto x,y,z, and reorder
    # nn will have structure
    #  n1x n1y n1z
    #  n2x n2y n2z
    #  n3x n3y n3z
    nn = 0.5 * abs(numpy.sqrt(-1.j * vn.conj().T.dot(jmt).dot(_vxyz[dms - 1])))
    rows = list(select)
    order = []
    for ixz in select:
        ind = numpy.argmax(nn[rows, ixz])
        order.append(rows[ind])
        del rows[ind]
    v_ordered = vn[:, order]
    lmbd = lmbd[order]
    aa = numpy.vstack((numpy.real(v_ordered), numpy.imag(v_ordered))).reshape(
        (nv, nv), order='F')
    return aa, lmbd


def amat(tt):
    """A = amat(T)
    Find the A matrix from one turn map matrix T
    Provided for backward compatibility, see " A, eigval = a_matrix(T)"

    INPUT
    T       (m, m)  transfer matrix for 1 turn

    OUTPUT
    A       (m, m)  A-matrix
    """
    aa, _ = a_matrix(tt)
    return aa


# noinspection PyPep8Naming
def symplectify(M):
    """
    symplectify makes a matrix more symplectic
    follow Healy algorithm as described by McKay
    BNL-75461-2006-CP
    """
    nv = M.shape[0]
    J = jmat(nv // 2)

    V = numpy.dot(J.dot(numpy.identity(nv) - M), inv(numpy.identity(nv) + M))
    # V should be almost symmetric.  Replace with symmetrized version.

    W = (V + V.T) / 2
    # Now reconstruct M from W
    JW = numpy.dot(J, W)
    MS = numpy.dot(numpy.identity(nv) + JW, inv(numpy.identity(nv) - JW))
    return MS


def get_mode_matrices(a):
    """Given a (m, m) A matrix , find the R-matrices of the m/2 normal modes"""

    def mul2(slc):
        return a[:, slc] @ tt[slc, slc]

    dms = a.shape[0] // 2
    # Rk = A * Ik * A.T                     Only for symplectic
    # modelist = [numpy.dot(a[:, sl], a.T[sl, :]) for sl in _submat[:dms]]
    # Rk = A * S * Ik * inv(A) * S.T        Even for non-symplectic
    ss = jmat(dms)
    tt = jmatswap(dms)
    a_s = numpy.concatenate([mul2(slc) for slc in _submat[:dms]], axis=1)
    inva = solve(a, ss.T)
    modelist = [a_s[:, sl] @ inva[sl, :] for sl in _submat[:dms]]
    return numpy.stack(modelist, axis=0)


def get_tunes_damp(tt, rr=None):
    """
    mode_emit, damping_rates, tunes = get_tunes_damp(T, R)

    INPUT
        T                   (m, m) transfer matrix for 1 turn
        R                   (m, m) beam matrix (optional)

        m can be 2 (single plane), 4 (betatron motion) or 6 (full motion)

    OUTPUT
        record array with the following fields:
        tunes               (m/2,) tunes of the m/2 normal modes
        damping_rates       (m/2,) damping rates of the m/2 normal modes
        mode_matrices       (m/2, m, m) the R-matrices of the m/2 normal modes
        mode_emittances     Only if R is specified: (m/2,) emittance of each
                            of the m/2 normal modes
    """
    nv = tt.shape[0]
    dms = int(nv / 2)
    aa, vps = a_matrix(tt)
    tunes = numpy.mod(numpy.angle(vps) / 2.0 / pi, 1.0)
    damping_rates = -numpy.log(numpy.absolute(vps))

    if rr is None:
        return numpy.rec.fromarrays(
            (numpy.array(tunes), numpy.array(damping_rates),
             numpy.array(get_mode_matrices(aa))),
            dtype=[('tunes', numpy.float64, (dms,)),
                   ('damping_rates', numpy.float64, (dms,)),
                   ('mode_matrices', numpy.float64, (dms, nv, nv))]
        )
    else:
        jmt = jmat(dms)
        rdiag = numpy.diag(aa.T.dot(jmt.dot(rr.dot(jmt.dot(aa)))))
        mode_emit = -0.5 * (rdiag[0:nv:2] + rdiag[1:nv:2])
        return numpy.rec.fromarrays(
            (numpy.array(tunes), numpy.array(damping_rates),
             numpy.array(get_mode_matrices(aa)), mode_emit),
            dtype=[('tunes', numpy.float64, (dms,)),
                   ('damping_rates', numpy.float64, (dms,)),
                   ('mode_matrices', numpy.float64, (dms, nv, nv)),
                   ('mode_emittances', numpy.float64, (dms,))]
        )
