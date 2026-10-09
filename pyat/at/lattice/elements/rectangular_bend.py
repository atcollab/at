"""Additional method for rectangular bending magnets"""

from math import sin, cos

import numpy as np
from scipy.optimize import fsolve

from .magnet_elements import Dipole

__all__ = []


def rbendtune(self: Dipole) -> None:
    # noinspection PyUnresolvedReferences
    r"""Set *X0ref* and *RefDZ* for rectangular bending magnets

    This method must be called after creating a rectangular bending magnet
    or after setting its *PolynomA/B* attributes. It will set the correct *X0ref*
    and *RefDZ* attributes to get a zero closed orbit for the reference particle.

    The method will do nothing on dipoles with a non-rectangular passmethod.

    Example:

        >>> # Identify the rectangular bends
        >>> rbends = ring.get_bool_index(...)
        >>> # Set their correct attributes
        >>> for dip in ring.select(rbends):
        ...     dip.rbendtune()

    """

    def cross(x0r: float):
        """Return the horizontal exit angle of the reference particle"""
        elem.X0ref = x0r
        out = elem.track(np.zeros(6))
        return out[1]

    def checkmul(el):
        """Check if there are multipoles"""
        for order in range(el.MaxOrder + 1):
            if el.PolynomB[order] != 0.0:
                return True
        return False

    passmethod = self.PassMethod.replace("RadPass", "Pass").replace("QuantPass", "Pass")
    if passmethod in {
        "BndStrMPoleSymplectic4Pass",
        "ExactRectangularBendPass",
        "ExactRectBendPass",
    }:
        elem = self.copy()
        elem.PassMethod = passmethod
        theta = elem.BendingAngle

        # Analytical estimate
        x0ref = elem.Length * ((cos(0.5 * theta) - 1.0) / theta + sin(0.5 * theta) / 12)

        # Search if there are multipoles
        if checkmul(self):
            # NumPy >=2.0 no longer allows float() on a non-0-d array; fsolve
            # returns shape (1,) for a scalar problem, so extract the element.
            x0ref = float(fsolve(cross, x0ref)[0])

        self.X0ref = x0ref

        # RefDZ must reflect the actual (deterministic) path length, including
        # mean radiation loss where applicable. A stochastic Quant passmethod
        # would otherwise inject one random photon-emission draw into a
        # reference/design attribute; substitute its deterministic RadPass
        # counterpart (validated to reproduce the Quant ensemble mean).
        refelem = self if not self.PassMethod.endswith("QuantPass") else self.copy()
        if refelem is not self:
            refelem.PassMethod = self.PassMethod.replace("QuantPass", "RadPass")
        # Radiating passmethods need Energy passed explicitly to track() -- the
        # element's own Energy attribute is not read automatically.
        energy = getattr(refelem, "Energy", None)
        rout = refelem.track(np.zeros(6), **({} if energy is None else {"energy": energy}))
        self.RefDZ = rout[5]


Dipole.rbendtune = rbendtune
