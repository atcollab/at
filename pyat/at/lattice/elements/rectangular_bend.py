"""Additional method for rectangular bending magnets."""

from math import sin, cos

import numpy as np
from scipy.optimize import fsolve

from .magnet_elements import Dipole

__all__ = []


def rbendtune(self: Dipole) -> None:
    # noinspection PyUnresolvedReferences
    r"""Set *X0ref* and *RefDZ* for rectangular bending magnets.

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
        """Return the horizontal exit angle of the reference particle."""
        elem.X0ref = x0r
        out = elem.track(np.zeros(6))
        return out[1]

    passmethod = self.PassMethod.replace("RadPass", "Pass")
    if passmethod in {
        "BndStrMPoleSymplectic4Pass",
        "ExactRectangularBendPass",
        "ExactRectBendPass",
    }:
        # check if there are multipoles
        if any(self.PolynomB[order] != 0.0 for order in range(self.MaxOrder + 1)):
            elem = self.copy()
            elem.PassMethod = passmethod
            tta = elem.BendingAngle

            # Analytical estimate
            x0ref = elem.Length * ((cos(0.5 * tta) - 1.0) / tta + sin(0.5 * tta) / 12)

            # cancel output angle
            x0ref = float(fsolve(cross, x0ref))
        else:
            x0ref = 0.0

        self.X0ref = x0ref
        rout = self.track(np.zeros(6))
        self.RefDZ = rout[5]


Dipole.rbendtune = rbendtune
