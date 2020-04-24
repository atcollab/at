import collections
import logging as log
from os.path import abspath
import re
import numpy

from at.lattice.elements import (
    Drift,
    Dipole,
    Quadrupole,
    Sextupole,
    Multipole,
    Corrector,
    Marker,
    RFCavity,
)
from at.lattice import Lattice
from at.load import register_format


def create_drift(name, params, energy, harmonic_number):
    length = parse_float(params.pop("l", 0))
    return Drift(name, length, **params)


def create_marker(name, params, energy, harmonic_number):
    return Marker(name, **params)


def create_quad(name, params, energy, harmonic_number):
    length = parse_float(params.pop("l", 0))
    params["NumIntSteps"] = params.pop("n", 10)
    params["k"] = parse_float(params.pop("k"))
    params["PassMethod"] = "StrMPoleSymplectic4Pass"
    return Quadrupole(name, length, **params)


def create_sext(name, params, energy, harmonic_number):
    length = parse_float(params.pop("l", 0))
    params["NumIntSteps"] = params.pop("n", 10)
    k2 = parse_float(params.pop("k", 0))
    return Sextupole(name, length, k2, **params)


def create_dipole(name, params, energy, harmonic_number):
    length = parse_float(params.pop("l", 0))
    params["NumIntSteps"] = params.pop("n", 10)
    params["PassMethod"] = "BndMPoleSymplectic4Pass"
    params["BendingAngle"] = (float(params.pop("t")) / 180) * numpy.pi
    params["EntranceAngle"] = (float(params.pop("t1")) / 180) * numpy.pi
    params["ExitAngle"] = (float(params.pop("t2")) / 180) * numpy.pi
    # Tracy is encoding gap plus fringe int in the 'gap' field.
    # Since BndMPoleSymplectic4Pass only uses the product of FringeInt
    # and gap we can substitute the following.
    if "gap" in params:
        params["FullGap"] = float(params.pop("gap")) / 2
        params["FringeInt1"] = 1
        params["FringeInt2"] = 1
    k = parse_float(params.pop("k", 0))
    if "PolynomB" in params:
        params["PolynomB"][1] = k
    else:
        params["PolynomB"] = [0, k, 0, 0]
    return Dipole(name, length, **params)


def create_corrector(name, params, energy, harmonic_number):
    length = parse_float(params.pop("l", 0))
    # Is there a property for this?
    kick_angle = [0, 0]
    return Corrector(name, length, kick_angle, **params)


def create_multipole(name, params, energy, harmonic_number):
    length = params.pop("l", 0)
    params["NumIntSteps"] = params.pop("n", 10)
    poly_a = [0, 0, 0, 0]
    poly_b = [0, 0, 0, 0]
    return Multipole(name, length, poly_a, poly_b, **params)


def create_cavity(name, params, energy, harmonic_number):
    length = parse_float(params.pop("l", 0))
    voltage = parse_float(params.pop("voltage"))
    frequency = parse_float(params.pop("frequency"))
    params["Phi"] = parse_float(params.pop("phi"))
    return RFCavity(name, length, voltage, frequency, harmonic_number, energy, **params)


ELEMENT_MAP = {
    "drift": create_drift,
    "bending": create_dipole,
    "quadrupole": create_quad,
    "sextupole": create_sext,
    "multipole": create_multipole,
    "corrector": create_corrector,
    "marker": create_marker,
    "map": create_marker,
    "beampositionmonitor": create_marker,
    "cavity": create_cavity,
}


def tokenise_expression(expression):
    tokens = []
    current_token = ""
    for char in expression:
        if char in "+-*/":
            if not current_token or current_token and current_token[-1] == "e":
                current_token += char
            else:
                tokens.append(current_token)
                tokens.append(char)
                current_token = ""
        else:
            current_token += char

    if current_token:
        tokens.append(current_token)
    return tokens


def parse_float(expression):
    """Can of worms."""
    log.debug("parse_float {}".format(expression))
    try:
        return float(expression)
    except ValueError:
        tokens = tokenise_expression(expression)
        log.debug(tokens)
        if ("*") in tokens:
            return float(tokens[0]) * float(tokens[2])
        elif ("/") in tokens:
            return float(tokens[0]) / float(tokens[2])
        elif ("+") in tokens:
            return float(tokens[0]) + float(tokens[2])
        elif ("-") in tokens:
            return float(tokens[0]) - float(tokens[2])


def split_ignoring_parentheses(string, delimiter):
    PLACEHOLDER = "placeholder"
    substituted = string[:]
    matches = collections.deque(re.finditer("\\(.*?\\)", string))
    for match in matches:
        substituted = substituted.replace(match.group(), PLACEHOLDER, 1)
    parts = substituted.split(delimiter)
    replaced_parts = []
    for part in parts:
        if PLACEHOLDER in part:
            next_match = matches.popleft()
            part = part.replace(PLACEHOLDER, next_match.group(), 1)
        replaced_parts.append(part)
    assert not matches

    return replaced_parts


def parse_lines(contents):
    """Return individual lines.

    Remove comments and whitespace, convert to lowercase, and split on
    semicolons.
    """
    # Nested comments not handled.
    in_comment = False
    stripped_contents = ""
    for char in contents.lower():
        if char == "{":
            in_comment = True
        elif char == "}":
            in_comment = False
        elif not in_comment and not char.isspace():
            stripped_contents += char
    return stripped_contents.split(";")


def parse_chunk(value, elements, chunks):
    log.debug('parse_chunk %s', value)
    chunk = []
    parts = split_ignoring_parentheses(value, ",")
    log.debug(parts)
    for part in parts:
        if "symmetry" in part:
            continue
        if "inv" in part:
            chunk_to_invert = re.match("inv\\((.*)\\)", part).groups()[0]
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
            num, contents = part.split("*")
            if contents.startswith("("):
                assert contents[-1] == ")"
                contents = contents[1:-1]
            chunk.extend(int(num) * parse_chunk(contents, elements, chunks))
        elif part in elements:
            chunk.append(elements[part])
        elif part in chunks:
            chunk.extend(chunks[part])
        else:
            raise Exception("part {} not understood".format(part))
    return chunk


def expand_tracy(contents, lattice_key, harmonic_number):
    lines = parse_lines(contents)
    assert lines[0] == "definelattice"
    assert lines[-2] == "end"
    assert lines[-1] == ""
    variables = {"harmonic_number": harmonic_number}
    elements = {}
    chunks = {}
    for line in lines[1:-2]:
        if ":" not in line:
            key, value = line.split("=")
            if key == "energy":
                value = parse_float(value) * 1e9
            variables[key] = value
        else:
            key, value = line.split(":")
            if value.split(",")[0].strip() in ELEMENT_MAP:
                elements[key] = tracy_element_from_string(key, value, variables)
            else:
                chunk = parse_chunk(value, elements, chunks)
                chunks[key] = chunk

    return chunks[lattice_key], variables["energy"]


def parse_hom(hom_string):
    """Parse 'hom' string from lattice file. Note that PolynomB
    is before PolynomA.

    hom(3, 1.2, 3.4) => [0, 0, 3.4], [0, 0, 1.2]
    """
    hom_parts = [part.strip() for part in hom_string.split(",")]
    assert len(hom_parts) % 3 == 0
    max_order = max(int(p) for i, p in enumerate(hom_parts) if i % 3 == 0)
    polynom_a = [0] * int(max_order)
    polynom_b = [0] * int(max_order)
    for i in range(len(hom_parts) // 3):
        order = int(hom_parts[i * 3])
        # PolynomB before PolynomA.
        b = parse_float(hom_parts[i * 3 + 1])
        a = parse_float(hom_parts[i * 3 + 2])
        polynom_a[order - 1] = a
        polynom_b[order - 1] = b

    return polynom_a, polynom_b


def tracy_element_from_string(name, element_string, variables):
    log.debug("Parsing tracy element {}".format(element_string))
    parts = split_ignoring_parentheses(element_string, ",")
    params = {}
    element_type = parts[0]
    for part in parts[1:]:
        try:
            key, value = split_ignoring_parentheses(part, "=")
        except ValueError:
            key, value = part, None
        if value in variables:
            value = variables[value]
        if key == "hom":
            assert value[0] == "("
            assert value[-1] == ")"
            polynom_a, polynom_b = parse_hom(value[1:-1])
            params["PolynomA"] = polynom_a
            params["PolynomB"] = polynom_b
        else:
            params[key] = value

    energy = variables["energy"]
    harmonic_number = variables["harmonic_number"]

    return ELEMENT_MAP[element_type](name, params, energy, harmonic_number)


def load_tracy(filename, **kwargs):
    try:
        harmonic_number = kwargs.pop("harmonic_number")
        lattice_key = kwargs.pop("lattice_key", "cell")

        def elem_iterator(params, tracy_file):
            with open(params.setdefault("tracy_file", tracy_file)) as f:
                contents = f.read()
                element_lines, energy = expand_tracy(
                    contents, lattice_key, harmonic_number
                )
                params.setdefault("energy", energy)
                for line in element_lines:
                    yield line

        return Lattice(abspath(filename), iterator=elem_iterator, **kwargs)
    except Exception as e:
        raise ValueError("Failed to load tracy lattice {}: {}".format(filename, e))


register_format(".lat", load_tracy, descr="Tracy format")
