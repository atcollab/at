"""Load a lattice from an Elegant file (.lte).

This is not complete but can parse the example files that I have.

Note that Elegant scales magnet polynomials in a different way
to AT, so the parsed coefficients need to be divided by n! for
the coefficient of order n.

"""
import collections
import logging as log
from os.path import abspath
import re

from at.lattice.elements import (
    Corrector,
    Drift,
    Dipole,
    Marker,
    Octupole,
    Quadrupole,
    RFCavity,
    Sextupole,
)
from at.lattice import Lattice
from at.load import register_format, utils


def create_drift(name, params, energy, harmonic_number):
    length = params.pop("l", 0)
    return Drift(name, length, **params)


def create_marker(name, params, energy, harmonic_number):
    return Marker(name, **params)


def create_quad(name, params, energy, harmonic_number):
    length = params.pop("l", 0)
    params["NumIntSteps"] = params.pop("n_kicks", 10)
    params["k"] = float(params.pop("k1"))
    params["PassMethod"] = "StrMPoleSymplectic4Pass"
    return Quadrupole(name, length, **params)


def create_sext(name, params, energy, harmonic_number):
    length = params.pop("l", 0)
    params["NumIntSteps"] = params.pop("n_kicks", 10)
    k2 = float(params.pop("k2", 0))
    return Sextupole(name, length, k2 / 2, **params)


def create_oct(name, params, energy, harmonic_number):
    length = params.pop("l", 0)
    params["NumIntSteps"] = params.pop("n_kicks", 10)
    k3 = float(params.pop("k3", 0))
    PolynomA = [0, 0, 0, 0]
    PolynomB = [0, 0, 0, k3]
    return Octupole(name, length, PolynomA, PolynomB, **params)


def create_dipole(name, params, energy, harmonic_number):
    length = params.pop("l", 0)
    params["NumIntSteps"] = params.pop("n_kicks", 10)
    params["PassMethod"] = "BndMPoleSymplectic4Pass"
    params["BendingAngle"] = float(params.pop("angle"))
    params["EntranceAngle"] = float(params.pop("e1"))
    params["ExitAngle"] = float(params.pop("e2"))
    if "hgap" in params:
        params["FullGap"] = float(params.pop("hgap")) * 2
        fint = params.pop("fint")
        params["FringeInt1"] = fint
        params["FringeInt2"] = fint
    k1 = float(params.pop("k1", 0))
    k2 = float(params.pop("k2", 0))
    k3 = float(params.pop("k3", 0))
    k4 = float(params.pop("k4", 0))
    params["PolynomB"] = [0, k1, k2 / 2, k3 / 6, k4 / 24]
    return Dipole(name, length, **params)


def create_corrector(name, params, energy, harmonic_number):
    length = params.pop("l", 0)
    hkick = params.pop("hkick", 0)
    vkick = params.pop("vkick", 0)
    kick_angle = [hkick, vkick]
    return Corrector(name, length, kick_angle, **params)


def create_cavity(name, params, energy, harmonic_number):
    length = params.pop("l", 0)
    voltage = params.pop("volt")
    frequency = params.pop("freq")
    params["Phi"] = params.pop("phase")
    return RFCavity(name, length, voltage, frequency, harmonic_number, energy, **params)


ELEMENT_MAP = {
    "drift": create_drift,
    "drif": create_drift,
    "edrift": create_drift,
    "csben": create_dipole,
    "csbend": create_dipole,
    "csrcsben": create_dipole,
    "quadrupole": create_quad,
    "kquad": create_quad,
    "ksext": create_sext,
    "koct": create_oct,
    "kicker": create_corrector,
    "hkick": create_corrector,
    "vkick": create_corrector,
    "mark": create_marker,
    "malign": create_marker,
    "recirc": create_marker,
    "sreffects": create_marker,
    "rcol": create_marker,
    "watch": create_marker,
    "charge": create_marker,
    "monitor": create_marker,
    "moni": create_marker,
    "rfca": create_cavity,
}


def parse_lines(contents):
    """Return individual lines.

    Remove comments and whitespace and convert to lowercase.
    """
    lines = []
    current_line = ""
    for line in contents.splitlines():
        line = line.lower()
        line = line.split("!")[0]
        line = line.strip()
        if line:
            if line.endswith("&"):
                current_line += line[:-1]
            else:
                lines.append(current_line + line)
                current_line = ""

    return lines


def parse_chunk(value, elements, chunks):
    """Parse a non-element part of a lattice file.

    line(a,b,c) => [a, b, c]
    """
    chunk = []
    parts = utils.split_ignoring_parentheses(value, ",")
    for part in parts:
        # What's this?
        if "symmetry" in part:
            continue
        if "line" in part:
            line_parts = re.match("line=\\((.*)\\)", part).groups()[0]
            for p in line_parts.split(","):
                p = p.strip()
                chunk.extend(parse_chunk(p, elements, chunks))
        elif part.startswith("-"):
            # Reverse a sequence of elements. When doing this, swap
            # the entrance and exit angles for any dipoles.
            chunk_to_invert = part[1:]
            inverted_chunk = []
            for el in reversed(chunks[chunk_to_invert]):
                if el.__class__ == Dipole:
                    inverted_dipole = el.copy()
                    inverted_dipole.EntranceAngle = el.ExitAngle
                    inverted_dipole.ExitAngle = el.EntranceAngle
                    inverted_chunk.append(inverted_dipole)
                else:
                    inverted_chunk.append(el)
            chunk.extend(inverted_chunk)
        elif "*" in part:
            num, chunk_name = part.split("*")
            if chunk_name[0] == "(":
                assert chunk_name[-1] == ")"
                chunk_name = chunk_name[1:-1]
            chunk.extend(int(num) * parse_chunk(chunk_name, elements, chunks))
        elif part in elements:
            chunk.append(elements[part])
        elif part in chunks:
            chunk.extend(chunks[part])
        else:
            raise Exception("part {} not understood".format(part))
    return chunk


def expand_elegant(contents, lattice_key, energy, harmonic_number):
    lines = parse_lines(contents)
    variables = {"energy": energy, "harmonic_number": harmonic_number}
    elements = {}
    chunks = {}
    for line in lines:
        if ":" not in line:
            key, value = line.split("=")
            variables[key.strip()] = value.strip()
        else:
            key, value = line.split(":")
            key = key.strip()
            value = value.strip()
            if value.split(",")[0] in ELEMENT_MAP:
                elements[key] = elegant_element_from_string(key, value, variables)
            else:
                chunk = parse_chunk(value, elements, chunks)
                chunks[key] = chunk

    return chunks[lattice_key.lower()]


def handle_value(value):
    value = value.strip()
    if value.startswith('"'):
        assert value.endswith('"')
        value = value[1:-1]
        value = value.split()
        if len(value) > 1:
            # Handle basic arithmetic e.g. "0.04 2 /" -> 0.02
            assert len(value) == 3
            if value[2] == "/":
                value = float(value[0]) / float(value[1])
            elif value[2] == "*":
                value = float(value[0]) * float(value[1])
            elif value[2] == "-":
                value = float(value[0]) - float(value[1])
            elif value[2] == "+":
                value = float(value[0]) + float(value[1])
        else:
            value = value[0]
    return value


def elegant_element_from_string(name, element_string, variables):
    log.debug("Parsing elegant element {}".format(element_string))
    parts = utils.split_ignoring_parentheses(element_string, ",")
    params = {}
    element_type = parts[0]
    for part in parts[1:]:
        key, value = utils.split_ignoring_parentheses(part, "=")
        key = key.strip()
        value = handle_value(value)
        if value in variables:
            value = variables[value]

        params[key] = value

    energy = variables["energy"]
    harmonic_number = variables["harmonic_number"]

    return ELEMENT_MAP[element_type](name, params, energy, harmonic_number)


def load_elegant(filename, **kwargs):
    try:
        energy = kwargs.pop("energy")
        lattice_key = kwargs.pop("lattice_key")
        harmonic_number = kwargs.pop("harmonic_number")

        def elem_iterator(params, elegant_file):
            with open(params.setdefault("elegant_file", elegant_file)) as f:
                contents = f.read()
                element_lines = expand_elegant(
                    contents, lattice_key, energy, harmonic_number
                )
                params.setdefault("energy", energy)
                for line in element_lines:
                    yield line

        lat = Lattice(abspath(filename), iterator=elem_iterator, **kwargs)
        return lat
    except Exception as e:
        raise ValueError("Failed to load elegant lattice {}: {}".format(filename, e))


register_format(".lte", load_elegant, descr="Elegant format")
