from __future__ import annotations

__all__ = ["Exporter"]

import sys
from collections.abc import Sequence, Generator

from .file_input import ElementDescr
from ..lattice import Lattice, elements as elt


class Exporter:
    delimiter: str = ";"
    continuation: str = ""
    label_fmt = str.maketrans("*/+-", "..__")  # Not allowed in destination format
    bool_fmt = None
    use_line = True

    def __init__(self, ring: Lattice, **kwargs):
        def store_elem(store, elem: ElementDescr):
            """Store an element in the selected dictionary"""

            def name_gen(name: str, max: int = 10000):
                yield name
                for i in range(1, max):
                    yield ".".join((name, str(i)))

            def check(name: str) -> tuple[bool, bool]:
                for st in self.all_stores:
                    if name in st:
                        if st is store:
                            if elem == st[name]:
                                return True, False
                            else:
                                return False, False
                        else:
                            return False, False
                return True, True

            for nm in name_gen(elem.name.translate(self.label_fmt)):
                valid, okstore = check(nm)
                if valid:
                    if okstore:
                        store[nm] = elem
                    return nm

            raise NameError(f"Cannot store {elem.name}")

        def scan(ring: Lattice) -> Generator[tuple[str, ElementDescr], None, None]:
            """Run through the lattice and store the converted elements"""
            end = 0.0
            for atelem in ring:
                attyp = type(atelem)
                elms = self.generate_madelems(attyp, vars(atelem).copy())
                if isinstance(elms, ElementDescr):
                    elms = (elms,)
                for elem in elms:
                    store = self.store.get(attyp, self.anystore)
                    length = elem.get("L", 0.0)
                    loc = end + 0.5 * length
                    end += length
                    if store is not None:
                        name = store_elem(store, elem)
                        yield name, loc

        use_line = kwargs.pop("use_line", self.use_line)
        self.seqname = kwargs.pop("use", ring.name if ring.name else "RING")
        self.anystore = {}
        self.store = {
            elt.Drift: {} if use_line else None,
            elt.Quadrupole: {},
            elt.Sextupole: {},
            elt.Dipole: {},
        }
        self.all_stores = [st for st in self.store.values() if st is not None] + [
            self.anystore
        ]
        ElementDescr.bool_fmt = self.bool_fmt
        self.seq = list(scan(ring))
        self.length = ring.cell_length
        self.in_file = getattr(ring, "in_file", "<unknown>")
        self.in_use = getattr(ring, "use", "<unknown")
        self.energy = ring.energy
        self.particle = ring.particle
        self.is_6d = ring.is_6d
        # self.anystore[self.seqname] = None

    def generate_madelems(
        self, eltype: type[elt.Element], elemdict: dict
    ) -> ElementDescr | list[ElementDescr]:
        pass

    def print_beam(self, file):
        pass

    def print_elems(self, store: dict, file) -> None:
        print(file=file)
        for elname, el in store.items():
            if el is not None:
                print(f"{elname.ljust(10)}: {str(el)}{self.delimiter}", file=file)

    def print_sequence(self, file) -> None:
        line = f"\n{self.seqname.ljust(10)}: SEQUENCE, L={self.length}{self.delimiter}"
        print(f"\n{line}{self.delimiter}", file=file)
        for elname, at in self.seq:
            print(f"  {elname.ljust(10)}, AT={at}{self.delimiter}", file=file)
        print(f"ENDSEQUENCE{self.delimiter}", file=file)

    def print_line(self, file) -> None:
        def elnames(subseq: Sequence):
            return ", ".join(ename.ljust(8) for ename, at in subseq)

        nl = int((len(self.seq) - 1) / 10)
        print(f"\n{self.seqname.ljust(10)}: LINE=( {self.continuation}", file=file)
        for lnum in range(nl):
            print(
                f"  {elnames(self.seq[10*lnum:10*(lnum+1)])}, {self.continuation}",
                file=file,
            )
        print(f"  {elnames(self.seq[10*nl:])}){self.delimiter}\n", file=file)

    def export(self, filename: str | None = None) -> None:
        def do_export(file):
            print(
                f"! Converted by PyAT from in_file={self.in_file}, use={self.in_use!r}\n",
                file=file,
            )
            self.print_beam(file)
            for store in self.all_stores:
                if store is not None:
                    self.print_elems(store, file)

            if self.store[elt.Drift] is None:
                self.print_sequence(file)
            else:
                self.print_line(file)

        if filename is None:
            do_export(sys.stdout)
        else:
            with open(filename, "w") as mfile:
                do_export(mfile)
