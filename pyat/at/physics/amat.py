"""A-matrix construction."""

from math import pi

import numpy as np
from scipy.linalg import block_diag, eig, inv, solve

from at.lattice import AtError

__all__ = [
    "a_matrix",
    "amat",
    "get_mode_matrices",
    "get_tunes_damp",
    "jmat",
    "jmatswap",
    "symplectify",
]

_i2 = np.array([[-1.0j, -1.0], [1.0, 1.0j]])

# Prepare symplectic identity matrix
_j2 = np.array([[0.0, 1.0], [-1.0, 0.0]])
_jm = [_j2, block_diag(_j2, _j2), block_diag(_j2, _j2, _j2)]
_jmswap = [_j2, block_diag(_j2, _j2), block_diag(_j2, _j2, _j2.T)]

_vxyz = [_i2, block_diag(_i2, _i2), block_diag(_i2, _i2, _i2)]
_submat = [slice(0, 2), slice(2, 4), slice(4, 6)]


def jmat(ind: int):
    """antisymetric block diagonal matrix [[0, 1][-1, 0]].

    Parameters:
        ind:    Matrix dimension: 1, 2 or 3

    Returns:
        S:      Block diagonal matrix, (2, 2) or (4, 4) or (6, 6)
    """
    return _jm[ind - 1]


def jmatswap(ind: int):
    """Modified antisymetric block diagonal matrix to deal with the swap of the
    longitudinal coordinates.
    """
    return _jmswap[ind - 1]


# noinspection PyPep8Naming
def a_matrix(M):
    r"""Find the :math:`\mathbf{A}` matrix from one turn map :math:`\mathbf{M}`.

    :math:`\mathbf{A}` represents a change of referential which converts the
    one-turn transfer matrix :math:`\mathbf{M}` into a set of rotations:

    .. math:: \mathbf{M} = \mathbf{A} \cdot \mathbf{R} \cdot \mathbf{A}^{-1}

    with :math:`\mathbf{R}` A block-diagonal matrix:

    .. math::

        \mathbf{R}=\begin{pmatrix}rot_1 & 0 & 0 \\ 0 & rot_2 & 0 \\ 0 & 0 &
        rot_3\end{pmatrix} \text{, and } rot_i = \begin{pmatrix}\cos{\mu_i} &
        \sin{\mu_i} \\ -\sin{\mu_i} & cos{\mu_i}\end{pmatrix}

    With radiation damping, the diagonal blocks are instead damped rotations:

    .. math::

        rot_i = \begin{pmatrix}\exp{(-\alpha_i}) & 0 \\ 0 & \exp{(-\alpha_i)}
        \end{pmatrix} \cdot \begin{pmatrix}\cos{\mu_i} & \sin{\mu_i} \\
        -\sin{\mu_i} & cos{\mu_i}\end{pmatrix}

    The order of diagonal blocks it set so that it is close to the order of
    x,y,z.

    Parameters:
        M:     (m, m)  transfer matrix for 1 turn

    m, the dimension of :math:`\mathbf{M}`, may be 2 (single plane),
    4 (betatron motion) or 6 (full motion)


    Returns:
        A:      (m, m)  A-matrix
        eigval: (m/2,)    Vector of Eigen values of T

    References:
        **[1]** Etienne Forest, Phys. Rev. E 58, 2481 – Published 1 August 1998
    """
    nv = M.shape[0]
    dms = int(nv / 2)
    jmt = jmatswap(dms)
    select = np.arange(0, nv, 2)
    rbase = np.stack((select, select), axis=1).flatten()

    # noinspection PyTupleAssignmentBalance
    lmbd, vv = eig(M)
    # Compute the norms
    vp = vv.conj().T @ jmt
    n = -0.5j * np.sum(vp.T * vv, axis=0)
    if any(abs(n) < 1.0e-12):
        msg = "Unstable ring"
        raise AtError(msg)
    # Move positive before negatives
    order = rbase + (n < 0)
    vv = vv[:, order]
    n = n[order]
    lmbd = lmbd[order]
    # Normalize vectors
    vn = vv / np.sqrt(abs(n)).reshape((1, nv))
    # find the vectors that project most onto x,y,z, and reorder
    # nn will have structure
    #  n1x n1y n1z
    #  n2x n2y n2z
    #  n3x n3y n3z
    nn = 0.5 * abs(np.sqrt(-1.0j * vn.conj().T @ jmt @ _vxyz[dms - 1]))
    rows = list(select)
    order = []
    for ixz in select:
        ind = np.argmax(nn[rows, ixz])
        order.append(rows[ind])
        del rows[ind]
    v_ordered = vn[:, order]
    lmbd = lmbd[order]
    aa = np.vstack((np.real(v_ordered), np.imag(v_ordered))).reshape(
        (nv, nv), order="F"
    )
    return aa, lmbd


# noinspection PyPep8Naming
def amat(M):
    """Find the A matrix from one turn map matrix T.

    Provided for backward compatibility, see :py:func:`a_matrix`

    Parameters:
        M:     (m, m)  transfer matrix for 1 turn

    Returns:
        A:     (m, m)  A-matrix
    """
    aa, _ = a_matrix(M)
    return aa


# noinspection PyPep8Naming
def symplectify(M):
    """Makes A matrix more symplectic.

    following the Healy algorithm described by MacKay

    Parameters:
        M:  Almost symplectic matrix

    Returns:
        MS: Symplectic matrix

    References:
        **[1]** `W.W.MacKay, Comment on Healy's symplectification algorithm,
        Proceedings of EPAC 2006
        <https://accelconf.web.cern.ch/e06/PAPERS/WEPCH152.PDF>`_
    """
    nv = M.shape[0]
    S = jmat(nv // 2)
    I = np.identity(nv)  # noqa: E741

    V = S @ (I - M) @ inv(I + M)
    # V should be almost symmetric.  Replace with symmetrised version.
    W = (V + V.T) / 2
    # Now reconstruct M from W
    SW = S @ W
    MS = (I + SW) @ inv(I - SW)
    return MS


# noinspection PyPep8Naming
def get_mode_matrices(A):
    """Derives the R-matrices from the A-matrix.

    Parameters:
        A:  A-matrix

    Returns:
        R:

    References:
        **[1]** Andrzej Wolski, Phys. Rev. ST Accel. Beams 9, 024001 –
        Published 3 February 2006
    """

    def mul2(slc):
        return A[:, slc] @ tt[slc, slc]

    dms = A.shape[0] // 2
    # Rk = A * Ik * A.T                     Only for symplectic
    # modelist = [np.dot(A[:, sl], A.T[sl, :]) for sl in _submat[:dms]]
    # Rk = A * S * Ik * inv(A) * S.T        Even for non-symplectic
    ss = jmat(dms)
    tt = jmatswap(dms)
    a_s = np.concatenate([mul2(slc) for slc in _submat[:dms]], axis=1)
    inva = solve(A, ss.T)
    modelist = [a_s[:, sl] @ inva[sl, :] for sl in _submat[:dms]]
    return np.stack(modelist, axis=0)


# noinspection PyPep8Naming
def get_tunes_damp(M, R=None):
    r"""Computes the mode emittances, tunes and damping times.

    Parameters:
        M:     (m, m) transfer matrix for 1 turn
        R:     (m, m) beam matrix (optional), allows computing the mode
                emittances

    m, the dimension of :math:`\mathbf{M}`, may be 2 (single plane),
    4 (betatron motion) or 6 (full motion)

    Returns:
        V:     record array with the following fields:

          **tunes**             (m/2,) tunes of the m/2 normal modes

          **damping_rates**     (m/2,) damping rates of the m/2 normal modes

          **mode_matrices**     (m/2, m, m) the R-matrices of the m/2 normal
          modes

          **mode_emittances**   Only if R is specified: (m/2,) emittance of each
          of the m/2 normal modes
    """
    nv = M.shape[0]
    dms = int(nv / 2)
    A, vps = a_matrix(M)
    tunes = np.mod(np.angle(vps) / 2.0 / pi, 1.0)
    damping_rates = -np.log(np.absolute(vps))

    if R is None:
        return np.rec.fromarrays(
            (
                np.array(tunes),
                np.array(damping_rates),
                np.array(get_mode_matrices(A)),
            ),
            dtype=[
                ("tunes", np.float64, (dms,)),
                ("damping_rates", np.float64, (dms,)),
                ("mode_matrices", np.float64, (dms, nv, nv)),
            ],
        )
    else:
        inva = inv(A)
        rdiag = np.diag(inva @ R @ inva.T)
        mode_emit = 0.5 * (rdiag[0:nv:2] + rdiag[1:nv:2])
        return np.rec.fromarrays(
            (
                np.array(tunes),
                np.array(damping_rates),
                np.array(get_mode_matrices(A)),
                mode_emit,
            ),
            dtype=[
                ("tunes", np.float64, (dms,)),
                ("damping_rates", np.float64, (dms,)),
                ("mode_matrices", np.float64, (dms, nv, nv)),
                ("mode_emittances", np.float64, (dms,)),
            ],
        )
