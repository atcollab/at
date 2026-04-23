"""Description of particles."""

import numpy as np

from ..constants import e_mass, p_mass, mu_mass


class Particle:
    """
    Particle object: it defines the properties of the particles circulating in a ring.
    """

    _known = {
        "relativistic": {"rest_energy": 0.0, "charge": -1.0},
        "electron": {"rest_energy": e_mass, "charge": -1.0},
        "positron": {"rest_energy": e_mass, "charge": 1.0},
        "proton": {"rest_energy": p_mass, "charge": 1.0},
        "antiproton": {"rest_energy": mu_mass, "charge": -1.0},
        "posmuon": {"rest_energy": mu_mass, "charge": 1.0},
        "negmuon": {"rest_energy": mu_mass, "charge": -1.0},
    }

    def __init__(self, name: str | None = "relativistic", **kwargs):
        """

        Parameters:
            name:   Particle name. 'electron', 'positron', 'proton', 'posmuon',
              'negmuon' are predefined. For other particles, the rest energy
              and charge must be provided as keywords. If *name* is blank or None, an
              attempt to identify the particle is done.

        Keyword Arguments:
            rest_energy:    Particle rest energy [ev]
            charge:         Particle charge [elementary charge]
            *:              Other keywords will be set as attributes of the
                              particle
        """

        def scan(dct):
            """Try to identify an unknown particle."""
            energy = dct["rest_energy"]
            charge = int(dct["charge"])
            for nm, val in self._known.items():
                de = abs(energy - val["rest_energy"]) / energy
                if de < 1.0e-8 and charge == int(val["charge"]):
                    return nm
            return "unknown"

        if name in self._known:
            kwargs.update(self._known[name])
        elif not name:
            # If name is blank, try to identify the particle
            name = scan(kwargs)

        self.name = name
        # Use a numpy scalar to allow division by zero
        self._rest_energy = np.array(kwargs.pop("rest_energy"), dtype=float)
        self._charge = kwargs.pop("charge")
        for key, value in kwargs.items():
            setattr(self, key, value)

    def to_dict(self) -> dict:
        attrs = vars(self).copy()
        attrs["rest_energy"] = attrs.pop("_rest_energy")
        attrs["charge"] = attrs.pop("_charge")
        return attrs

    def __repr__(self):
        if self.name in self._known:
            return f"Particle('{self.name}')"
        else:
            attrs = self.to_dict()
            name = attrs.pop("name")
            args = ", ".join(f"{k}={v!r}" for k, v in attrs.items())
            return f"Particle('{name}', {args})"

    def __str__(self):
        if self.name in self._known:
            return self.name
        else:
            return self.__repr__()

    # Use properties so that they are read-only
    @property
    def rest_energy(self) -> np.ndarray:
        """Particle rest energy [eV]."""
        return self._rest_energy

    @property
    def charge(self) -> float:
        """Particle charge [elementary charge]."""
        return self._charge
