"""Crab cavity Element"""

from __future__ import annotations

from .conversions import _array
from .abstract_elements import LongtMotion
from .basic_elements import LongElement


class CrabCavity(LongtMotion, LongElement):
    """Crab cavity element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + [
        "Voltages",
        "Frequency",
        "HarmNumber",
    ]
    default_pass = {False: "DriftPass", True: "CrabCavityPass"}
    _conversions = dict(
        LongElement._conversions,
        Voltages=lambda v: _array(v, (2,)),
        Frequency=float,
        HarmNumber=int,
        TimeLag=float,
        PhaseLag=float,
        SigPhi=float,
        SigVV=float,
    )

    def __init__(
        self,
        family_name: str,
        length: float,
        voltages: float,
        frequency: float,
        harmonic_number: int,
        **kwargs,
    ):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            voltages:       [Horizontal voltage, Vertical voltage] [V]
            frequency:      RF frequency [Hz]
            harmonic_number:

        Keyword Args:
            TimeLag=0:      Cavity time lag [m]
            SigPhi=0:       Phase noise [rad]
            SigVV=0:        Voltage noise [V]

        Default PassMethod: ``CrabCavityPass``
        """
        kwargs.setdefault("TimeLag", 0.0)
        kwargs.setdefault("PassMethod", self.default_pass[True])
        super().__init__(
            family_name,
            length,
            Voltages=voltages,
            Frequency=frequency,
            HarmNumber=harmonic_number,
            **kwargs,
        )

    def _part(self, fr, sumfr):
        pp = super()._part(fr, sumfr)
        pp.Voltages = fr * self.Voltages
        return pp

    def _get_longt_motion(self):
        return self.PassMethod.endswith("CavityPass")

    # noinspection PyShadowingNames
    def set_longt_motion(self, enable, new_pass=None, **kwargs):
        if new_pass == "auto":
            new_pass = (
                self.default_pass[True]
                if enable
                else ("IdentityPass" if self.Length == 0 else "DriftPass")
            )
        return super().set_longt_motion(enable, new_pass=new_pass, **kwargs)
