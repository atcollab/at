from .elements import Dipole
import numpy as np
from scipy.optimize import fsolve
from math import sin, cos

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
        """Check if there are multipole"""
        for order in range(el.MaxOrder + 1):
            if el.PolynomB[order] != 0.0:
                return True
        return False

    passmethod = self.PassMethod.replace("RadPass", "Pass")
    if any(
        [
            passmethod == pm
            for pm in (
                "BndStrMPoleSymplectic4Pass",
                "ExactRectangularBendPass",
                "ExactRectBendPass",
            )
        ]
    ):
        elem = self.copy()
        elem.PassMethod = passmethod
        theta = elem.BendingAngle

        # Analytical estimate
        x0ref = elem.Length * ((cos(0.5 * theta) - 1.0) / theta + sin(0.5 * theta) / 12)

        # Search if there are multipoles
        if checkmul(self):
            x0ref = float(fsolve(cross, x0ref))

        self.X0ref = x0ref
        rout = self.track(np.zeros(6))
        self.RefDZ = rout[5]


Dipole.rbendtune = rbendtune
