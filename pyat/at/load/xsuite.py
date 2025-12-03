from __future__ import annotations
from os.path import abspath
from typing import Optional, Any
import numpy as np
import copy
import warnings
from collections.abc import Sequence, Generator, Iterable
import json
from ..lattice import Particle, transformation
from ..lattice import Lattice, RFCavity
from ..lattice import Element, Drift, ThinMultipole, Marker
from ..lattice import AtWarning
from ..lattice import constants as cst

__all__ = ["save_xsuite", "load_xsuite"]

_multipole = ["ThinMultipole", "Multipole", "Quadrupole", "Sextupole", "Octupole"]

_dipole = ["Bend", "RBend", "Dipole"]

_cavity = ["Cavity", "RFCavity"]


class BasicElement(Element):

    _xsuite2at_attr = {
        "length": "Length",
        "__class__": "Class",
    }

    _xsuite2at_class = {
        "Quadrupole": "Quadrupole",
        "Sextupole": "Sextupole",
        "Multipole": "ThinMultipole",
        "RBend": "Dipole",
        "Bend": "Dipole",
        "Octupole": "Multipole",
        "Cavity": "RFCavity",
    }

    _at2xsuite_class = {
        "Quadrupole": "Quadrupole",
        "Sextupole": "Sextupole",
        "ThinMultipole": "Multipole",
        "Dipole": "Bend",
        "Bend": "Bend",
        "RFCavity": "Cavity",
    }

    def __init__(self, name: str, xsuite_params: dict = {}, atparams: dict = {}):
        self.name = name
        self.xsuite_params = copy.deepcopy(xsuite_params)
        self.atparams = copy.deepcopy(atparams)
        self._at2xsuite_attr = {v: k for k, v in self._xsuite2at_attr.items()}
        super().__init__(name)

    def _params_to_at(self):
        self.atparams["FamName"] = self.name
        if "length" not in self.xsuite_params.keys():
            self.xsuite_params["length"] = 0.0
        for k, v in self.xsuite_params.items():
            atk = self._xsuite2at_attr.get(k, None)
            if atk is not None:
                self.atparams[atk] = v

    def _params_to_xsuite(self):
        for k, v in self.atparams.items():
            xsk = self._at2xsuite_attr.get(k, None)
            if xsk is not None:
                self.xsuite_params[xsk] = self._convert_val(v)

    @staticmethod
    def _convert_val(v):
        if isinstance(v, np.ndarray):
            return v.tolist()
        elif isinstance(v, np.integer):
            return int(v)
        elif isinstance(v, np.floating):
            return float(v)
        else:
            return v

    @staticmethod
    def _get_default_class(length):
        defaultclass = "Marker"
        if length > 0.0:
            defaultclass = "Drift"
        return defaultclass

    def _class_to_at(self):
        length = self.atparams.get("Length", 0.0)
        defaultclass = self._get_default_class(length)
        cls = self.atparams["Class"]
        self.atparams["Class"] = self._xsuite2at_class.get(cls, defaultclass)

    def _class_to_xsuite(self):
        length = self.xsuite_params.get("length", 0.0)
        defaultclass = self._get_default_class(length)
        cls = self.xsuite_params["__class__"]
        self.xsuite_params["__class__"] = self._at2xsuite_class.get(cls, defaultclass)

    def get_at_element(self) -> Element:
        self._params_to_at()
        self._class_to_at()
        return super().from_dict(self.atparams)

    def get_xsuite_dict(self) -> dict:
        self._params_to_xsuite()
        self._class_to_xsuite()
        return self.xsuite_params


class MultipoleElement(BasicElement):

    _mag_order = {
        "Quadrupole": 1,
        "Sextupole": 2,
        "Octupole": 3,
    }

    _xsuite2at_attr = BasicElement._xsuite2at_attr
    _xsuite2at_attr.update(
        {
            "order": "MaxOrder",
            "num_multipole_kicks": "NumIntSteps",
        }
    )

    def __init__(self, name: str, xsuite_params: dict = {}, atparams: dict = {}):
        super().__init__(name, xsuite_params, atparams)

    def _set_at_transforms(self, elem):
        shift_x = self.xsuite_params.pop("shift_x", 0.0)
        shift_y = self.xsuite_params.pop("shift_y", 0.0)
        shift_s = self.xsuite_params.pop("shift_s", 0.0)
        rot_x = self.xsuite_params.pop("rot_x_rad", 0.0)
        rot_y = self.xsuite_params.pop("rot_y_rad", 0.0)
        rot_s = self.xsuite_params.pop("rot_s_rad_no_frame", 0.0)
        anchor = self.xsuite_params.pop("rot_shift_anchor", 0.0)
        if anchor == 0.0:
            reference = transformation.ReferencePoint.ENTRANCE
        elif anchor == 0.5:
            reference = transformation.ReferencePoint.CENTRE
        else:
            msg = "Anchor point for rotation different from 0 or 1, setting" \
            "to AT default: ReferencePoint.CENTRE"
            warnings.warn(AtWarning(msg))
            reference = transformation.ReferencePoint.CENTRE              
        
        if np.any(np.array([shift_x, shift_y, shift_s, rot_x, rot_y, rot_s]) != 0.0):
            transforms = {
                "dx": shift_x,
                "dy": shift_y,
                "dz": shift_s,
                "tilt": rot_s,
                "pitch": -rot_x,
                "yaw": rot_y,
                "reference": reference,
            }
            print(transforms)
            transformation.transform_elem(elem, **transforms)

    def _set_xs_transforms(self, dict_elem):
        transforms = self.atparams.get("transforms", None)
        if transforms is not None:
            dict_elem.update(transforms)

    def _set_at_poly(self):
        l = self.xsuite_params.get("length", 0.0)
        k1 = self.xsuite_params.pop("k1", 0.0)
        k1s = self.xsuite_params.pop("k1s", 0.0)
        k2 = self.xsuite_params.pop("k2", 0.0)
        k2s = self.xsuite_params.pop("k2s", 0.0)
        k3 = self.xsuite_params.pop("k3", 0.0)
        k3s = self.xsuite_params.pop("k3s", 0.0)
        knl = np.array(self.xsuite_params.pop("knl", np.zeros(4)))
        if len(knl) < 4:
            knl = np.concatenate((knl, np.zeros(4 - len(knl))))
        ksl = np.array(self.xsuite_params.pop("ksl", np.zeros(4)))
        if len(ksl) < 4:
            ksl = np.concatenate((ksl, np.zeros(4 - len(ksl))))
        if l > 0.0:
            knl /= l
            ksl /= l
        kn = np.concatenate(([0.0, k1, k2, k3], np.zeros(len(knl) - 4)))
        ks = np.concatenate(([0.0, k1s, k2s, k3s], np.zeros(len(ksl) - 4)))
        kn += knl
        ks += ksl
        self.atparams["PolynomB"] = self.p_to_at(list(kn))
        self.atparams["PolynomA"] = self.p_to_at(list(ks))
        if self.xsuite_params.get('angle', 0.0) == 0.0:
            if self.xsuite_params.get("edge_entry_active", 0) ==1:
                self.atparams["FringeQuadEntrance"] = 1
            if self.xsuite_params.get("edge_exit_active", 0) == 1:
                self.atparams["FringeQuadExit"] = 1

    def _set_xs_poly(self):
        pola = self.p_to_xsuite(self.atparams.pop("PolynomA", np.zeros(4)))
        polb = self.p_to_xsuite(self.atparams.pop("PolynomB", np.zeros(4)))
        l = self.atparams.get("Length")
        korder = self._mag_order.get(self.atparams.get("Class"), None)
        if korder is not None:
            self.xsuite_params["k" + str(korder)] = polb[korder]
            self.xsuite_params["k" + str(korder) + "s"] = pola[korder]
            pola[korder] = 0.0
            polb[korder] = 0.0
        if l > 0:
            polb *= l
            pola *= l
        if np.any(pola) or np.any(polb):
            self.xsuite_params["knl"] = list(polb)
            self.xsuite_params["ksl"] = list(pola)
        if korder == 1:
            qfrin = self.atparams.get("FringeQuadEntrance", 0)
            qfrout = self.atparams.get("FringeQuadExit", 0)
            self.xsuite_params["edge_entry_active"] = qfrin
            self.xsuite_params["edge_exit_active"] = qfrout              
        elif self.atparams.get('BendingAngle', 0.0) == 0.0:
            self.xsuite_params["edge_exit_active"] = 0
            self.xsuite_params["edge_entry_active"] = 0



    @staticmethod
    def poly_from_xsuite(x: Iterable[float], factor: float = 1.0) -> Generator[float]:
        """Convert polynomials from XSUITE to AT"""
        f = 1.0
        for n, vx in enumerate(x):
            yield factor * float(vx / f)
            f *= n + 1

    @staticmethod
    def poly_to_xsuite(x: Iterable[float], factor: float = 1.0) -> Generator[float]:
        """Convert polynomials from AT to XSUITE"""
        f = 1.0
        for n, vx in enumerate(x):
            yield factor * float(vx * f)
            f *= n + 1

    def p_to_at(self, a: float | Sequence[float]) -> np.ndarray:
        """Convert polynomials from XSUITE to AT"""
        if not isinstance(a, Sequence):
            # In case of a single element, we have a scalar instead of a tuple
            a = (a,)
        return np.fromiter(self.poly_from_xsuite(a), dtype=float)

    def p_to_xsuite(self, a: float | Sequence[float]) -> np.ndarray:
        """Convert polynomials from XSUITE to AT"""
        # if not isinstance(a, Sequence):
        #    # In case of a single element, we have a scalar instead of a tuple
        #    a = (a,)
        return np.fromiter(self.poly_to_xsuite(a), dtype=float)

    def get_at_element(self) -> Element:
        self._set_at_poly()
        elem = super().get_at_element()
        self._set_at_transforms(elem)
        return elem

    def get_xsuite_dict(self) -> dict:
        self._set_xs_poly()
        dict_elem = super().get_xsuite_dict()
        self._set_xs_transforms(dict_elem)
        return dict_elem


class DipoleElement(MultipoleElement):

    _xsuite2at_attr = MultipoleElement._xsuite2at_attr
    _xsuite2at_attr.update(
        {
            "angle": "BendingAngle",
            "edge_entry_fint": "FringeInt1",
            "edge_exit_fint": "FringeInt2",
            "edge_entry_angle": "EntranceAngle",
            "edge_exit_angle": "ExitAngle",
        }
    )

    def __init__(self, name: str, xsuite_params: dict = {}, atparams: dict = {}):
        super().__init__(name, xsuite_params, atparams)
        self.xsuite_params["model"] = "drift-kick-drift-exact"

    def _set_at_angle_faces(self):
        entry_hgap = self.xsuite_params.pop("edge_entry_hgap", None)
        exit_hgap = self.xsuite_params.pop("edge_exit_hgap", None)
        if entry_hgap != exit_hgap:
            msg = "Entry and Exit gaps for dipole are different, use entry"
            warnings.warn(AtWarning(msg))
        if entry_hgap is not None:
            self.atparams.update({"FullGap": entry_hgap})
        if self.xsuite_params.get("__class__") == "RBend":
            hangle = abs(0.5 * self.xsuite_params.get("angle"))
            self.xsuite_params["edge_entry_angle"] += hangle
            self.xsuite_params["edge_exit_angle"] += hangle
            self.xsuite_params["length"] /= np.sinc(hangle)

        entry_active = self.xsuite_params.get("_edge_entry_active", 1)
        entry_model = self.xsuite_params.get("_edge_entry_model", 0)
        if entry_model == -1 or entry_active==0:
            self.atparams["FringeBendEntrance"] = 0
        elif entry_model == 1:
            self.atparams["FringeQuadEntrance"] = 1
 
        exit_active = self.xsuite_params.get("_edge_exit_active", 1)
        exit_model = self.xsuite_params.get("_edge_exit_model", 0)
        if exit_model == -1 or exit_active==0:
            self.atparams["FringeBendExit"] = 0
        elif exit_model == 1:
            self.atparams["FringeQuadExit"] = 1

    def _set_xs_angle_faces(self):
        gap = self.atparams.pop("FullGap", None)
        if gap is not None:
            self.xsuite_params.update(
                {
                    "edge_entry_gap": gap,
                    "edge_exit_gap": gap,
                }
            )
        self.xsuite_params.update({"k0_from_h": True})
        dfr_entry = self.atparams.pop("FringeBendEntrance", 1)
        dfr_exit = self.atparams.pop("FringeBendExit", 1)
        qfr_entry = self.atparams.pop("FringeQuadEntrance", 0)
        qfr_exit = self.atparams.pop("FringeQuadExit", 0)

        if dfr_entry == 0:
            self.xsuite_params["edge_entry_model"] = "suppressed"
        elif qfr_entry == 1:
            self.xsuite_params["edge_entry_model"] = "full"
        else:
            self.xsuite_params["edge_entry_model"] = "linear"

        if dfr_exit == 0:
            self.xsuite_params["edge_exit_model"] = "suppressed"
        elif qfr_exit == 1:
            self.xsuite_params["edge_exit_model"] = "full"
        else:
            self.xsuite_params["edge_exit_model"] = "linear"

    def get_at_element(self) -> Element:
        self._set_at_angle_faces()
        return super().get_at_element()

    def get_xsuite_dict(self) -> dict:
        self._set_xs_angle_faces()
        return super().get_xsuite_dict()


class CavityElement(BasicElement):

    _xsuite2at_attr = BasicElement._xsuite2at_attr
    _xsuite2at_attr.update(
        {
            "voltage": "Voltage",
            "frequency": "Frequency",
        }
    )

    def __init__(
        self, name: str, xsuite_params: dict = {}, atparams: dict = {}, **kwargs
    ):
        super().__init__(name, xsuite_params, atparams)
        energy = kwargs.get("energy", 1.0e9)
        harmnumber = kwargs.get("harmnumber", 1)
        self.atparams["energy"] = energy
        self.atparams["harmonic_number"] = harmnumber

    def _set_at_lag(self):
        om = 2 * np.pi * self.xsuite_params.get("frequency")
        pl = self.xsuite_params.pop("lag", 0.0)
        tl = pl / om * cst.clight / 360
        self.atparams["TimeLag"] = tl

    def _set_xs_lag(self):
        om = 2 * np.pi * self.atparams.get("Frequency")
        tl = self.atparams.pop("TimeLag", 0.0)
        pl = tl * om / cst.clight * 360
        self.xsuite_params["lag"] = pl

    def get_at_element(self) -> Element:
        self._set_at_lag()
        return super().get_at_element()

    def get_xsuite_dict(self) -> dict:
        self._set_xs_lag()
        return super().get_xsuite_dict()


def at_from_xsuite(name: str, xsuite_params: dict = {}) -> BasicElement:
    cls = xsuite_params["__class__"]
    if cls in _dipole:
        elem = DipoleElement(name, xsuite_params=xsuite_params)
    elif cls in _multipole:
        elem = MultipoleElement(name, xsuite_params=xsuite_params)
    elif cls in _cavity:
        elem = CavityElement(name, xsuite_params=xsuite_params)
    else:
        elem = BasicElement(name, xsuite_params=xsuite_params)
    return elem.get_at_element()


def at_to_xsuite(name: str, atparams: dict = {}) -> dict:
    cls = atparams["Class"]
    if cls in _dipole:
        elem = DipoleElement(name, atparams=atparams)
    elif cls in _multipole:
        elem = MultipoleElement(name, atparams=atparams)
    elif cls in _cavity:
        elem = CavityElement(name, atparams=atparams)
    else:
        elem = BasicElement(name, atparams=atparams)
    return elem.get_xsuite_dict()


def load_xsuite(filename: str, **kwargs) -> Lattice:

    def json_generator(params: dict[str, Any], fn):
        with open(params.setdefault("in_file", fn), "rt") as jsonfile:
            data = json.load(jsonfile)

        elements = data["elements"]
        try:
            particle = data["particle_ref"]
            atparticle = Particle("relativistic")
            mass0 = float(particle["mass0"])
            gamma0 = float(particle["gamma0"][0])
            energy = float(gamma0 * mass0)
        except KeyError:
            atparticle = Particle("relativistic")
            energy = 1e9

        beam_params = {
            "particle": atparticle,
            "energy": energy,  # [eV]
            "periodicity": 1,
        }

        for k, v in beam_params.items():
            params.setdefault(k, v)
        for k, v in elements.items():
            yield at_from_xsuite(k, v)

    lattice = Lattice(abspath(filename), iterator=json_generator, **kwargs)

    hmin = np.inf
    for e in lattice[RFCavity]:
        e.HarmNumber = int(e.Frequency / lattice.revolution_frequency)
        if e.HarmNumber < hmin:
            hmin = e.HarmNumber
    lattice.harmonic_number = hmin

    return lattice


def _generate_thin(thick_multipole):
    d_thick = thick_multipole.to_dict()
    _ = d_thick.pop("PassMethod")
    name = d_thick.pop("FamName")
    length = d_thick.pop("Length")
    dr_entrance = Drift(name + "_Entrance", length / 2)
    dr_exit = Drift(name + "_Exit", length / 2)
    mult = ThinMultipole(
        name,
        d_thick.pop("PolynomA") * length,
        d_thick.pop("PolynomB") * length,
        **d_thick,
    )
    return [dr_entrance, mult, dr_exit]


def _set_unique_names(lattice):
    names = [e.FamName for e in lattice]
    _, idx_inv, cnt = np.unique(names, return_inverse=True, return_counts=True)
    idx_list = np.split(np.argsort(idx_inv), np.cumsum(cnt[:-1]))
    for ia in idx_list:
        if len(ia) > 1:
            for i, ii in enumerate(ia[::-1]):
                names[ii] += "_" + str(i)
    for e, name in zip(lattice, names):
        e.FamName = name


def _format_lattice(lattice):
    idx_mult = [i for i, e in enumerate(lattice) if e.__class__.__name__ == "Multipole"]
    if len(idx_mult) > 0:
        newring = lattice[: idx_mult[0]]
        for i, ii in enumerate(idx_mult[:-1]):
            newring += _generate_thin(lattice[ii])
            newring += lattice[ii + 1 : idx_mult[i + 1]]
        newring += _generate_thin(lattice[idx_mult[-1]])
        newring += lattice[idx_mult[-1] + 1 :]
    else:
        newring = lattice.deepcopy()
    _set_unique_names(newring)
    return Lattice(newring, energy=lattice.energy, periodicity=lattice.periodicity)


def save_xsuite(
    lattice: Lattice, filename: str, compact: bool = False, **kwargs
) -> None:
    lattice = _format_lattice(lattice)
    indent = None if compact else 2
    for e in lattice:
        refpoint = transformation.ReferencePoint[getattr(e, 'ReferencePoint', 'CENTRE')]
        offset, tilt, yaw, pitch = transformation.get_offsets_rotations(e, reference=refpoint)
        if np.any(np.array([offset[0], offset[1], offset[2], tilt, pitch, yaw]) != 0.0):
            if refpoint == transformation.ReferencePoint.ENTRANCE:
                anchor = 0
            else:
                anchor = 0.5
            e.transforms = {
                "shift_x": offset[0],
                "shift_y": offset[1],
                "shift_s": offset[2],
                "rot_s_rad_no_frame": tilt,
                "rot_x_rad": -pitch,
                "rot_y_rad": yaw,
                "rot_shift_anchor": anchor,
            }
    elements = {e.FamName: at_to_xsuite(e.FamName, e.to_dict()) for e in lattice}
    element_names = [e.FamName for e in lattice]
    with open(filename, "wt") as jsonfile:
        json.dump(
            {
                "elements": elements,
                "element_names": element_names,
            },
            jsonfile,
            indent=indent,
        )
