from __future__ import annotations

__all__ = ["translate"]

from collections.abc import Iterable, Sequence, Callable
import itertools

from .file_input import ElementDescr
from ..lattice import Lattice, elements as elt


class Ignore(ElementDescr):
    @classmethod
    def from_at(cls, kwargs):
        print(f"Ignore {kwargs['name']}")
        return []


class Beam(ElementDescr):
    trans_table = ()

    def __init__(self, ring: Lattice):
        super().__init__({"PassMethod": None, "FamName": ring.name})
        self["ENERGY"] = 1.0e-9 * ring.energy
        part = str(ring.particle)
        if part == "relativistic":
            part = "electron"
        self["PARTICLE"] = part.upper()
        self["RADIATE"] = ring.is_6d


def generate_madelems(
    trans: Callable[[type[elt.Element], type[ElementDescr]]],
    atelems: Iterable[elt.Element],
):
    for atel in atelems:
        attyp = type(atel)
        madel = trans(attyp)(vars(atel).copy())
        if isinstance(madel, ElementDescr):
            yield attyp, madel
        elif isinstance(madel, Iterable):
            yield from zip(itertools.repeat(attyp), madel)


def translate(
    trans: Callable[[type[elt.Element], type[ElementDescr]]],
    ring: Lattice,
    file,
    delimiter: str = ";",
    continuation: str = "",
    bool_fmt=None,
    use_line: bool = False,
    beam_descr: Callable | None = None,
):
    def store_elem(store, elem: ElementDescr):
        """Store an element in the selected dictionary"""

        def name_gen(name: str, max: int = 10000):
            yield name
            for i in range(1, max):
                yield ".".join((name, str(i)))

        def check(name: str) -> tuple[bool, bool]:
            for st in all_stores:
                if name in st:
                    if st is store:
                        if elem == st[name]:
                            return True, False
                        else:
                            return False, False
                    else:
                        return False, False
            return True, True

        for nm in name_gen(elem.name.translate(elem.label_fmt)):
            valid, okstore = check(nm)
            if valid:
                if okstore:
                    store[nm] = elem
                return nm

        raise NameError(f"Cannot store {elem.name}")

    def scan(ring, params):
        """Run through the lattice and store the converted elements"""
        end = 0.0
        for attyp, elem in generate_madelems(trans, ring):
            store = getstore.get(attyp, anystore)
            length = elem.get("L", 0.0)
            loc = end + 0.5 * length
            end += length
            if store is not None:
                name = store_elem(store, elem)
                yield name, loc
        params["L"] = end

    def print_elems(store: dict, file) -> None:
        print(file=file)
        for elname, el in store.items():
            if el is not None:
                print(f"{elname.ljust(10)}: {str(el)}{delimiter}", file=file)

    def print_sequence(seq: Sequence, seqname: str, length: float, file) -> None:
        print(f"\n{seqname.ljust(10)}: SEQUENCE, L={length}{delimiter}", file=file)
        for elname, at in seq:
            print(f"  {elname.ljust(10)}, AT={at}{delimiter}", file=file)
        print(f"ENDSEQUENCE{delimiter}", file=file)

    def print_line(seq: Sequence, seqname, file) -> None:
        def elnames(subseq: Sequence):
            return ", ".join(ename.ljust(8) for ename, at in subseq)

        nl = int((len(seq) - 1) / 10)
        print(f"\n{seqname.ljust(10)}: LINE=( {continuation}", file=file)
        for lnumber in range(nl):
            print(
                f"  {elnames(seq[10*lnumber:10*(lnumber+1)])}, {continuation}",
                file=file,
            )
        print(f"  {elnames(seq[10*nl:])}){delimiter}\n", file=file)

    ElementDescr.bool_fmt = bool_fmt
    anystore = {}
    getstore = {
        elt.Drift: {} if use_line else None,
        elt.Quadrupole: {},
        elt.Sextupole: {},
        elt.Dipole: {},
        Ignore: None,
    }
    all_stores = [st for st in getstore.values() if st is not None] + [anystore]

    seqname = ring.name if ring.name else "RING"
    anystore[seqname] = None
    params = {}
    seq = list(scan(ring, params))

    in_file = getattr(ring, "in_file", "<unknown>")
    in_use = getattr(ring, "use", "<unknown")
    print(f"! Converted by PyAT from in_file={in_file}, use={in_use!r}\n", file=file)
    if beam_descr:
        print(f"{beam_descr(ring)}{delimiter}", file=file)

    for dd in getstore.values():
        if dd is not None:
            print_elems(dd, file)
    print_elems(anystore, file)

    if use_line:
        print_line(seq, seqname, file)
    else:
        print_sequence(seq, seqname, params["L"], file)
