from math import sqrt  # noqa: F401

import numpy as np
import pytest
from numpy.testing import assert_allclose as assert_close

from at import UnorderedParser, MadxParser, Mad8Parser, ElegantParser, set_argparser

__all__ = []

test_data1 = """
/* test block comment
delimiter:      ";"
linecomment:    "!", "//"
end of comments */

v1 = 42.0 ;                  // Test line comment
v2 = 2*(sqrt(4) + 2**2) ;    // Test expression

label : command1, flag1,     ! Test continuation
-FLAG2, title="test parser", ! Test disabled flag
arg1=v1, ARG2=v2, ARG3=V3 ;  ! Test postponed definition

V3 = True ; V4 = False ;     ! Test several commands
"""

test_data2 = r"""
/* test block comment
delimiter:      end-of-line
linecomment:    "#"
continuation:   "\"
end of comments */

V1 = 42.0                     # Test line comment
V2 = 2*(sqrt(4) + 2**2)       # Test expression

label : command1, flag1, \    # Test continuation
-FLAG2, title="test parser",\ # Test disabled flag
arg1=V1, ARG2=V2, ARG3=V3     # Test postponed definition

V3 = True ; V4 = False        # Test several commands

TEST: TWISS, D=BETA;
"""

madx_data = """
/* test block comment
delimiter:      ";"
linecomment:    "!", "//"
end of comments */

BEAM,   PARTICLE="electron", ENERGY=6.0, SEQUENCE=RING.1;
BEAM,   PARTICLE="positron", ENERGY=2.0;
BEAM,   SEQUENCE=RING.1, RADIATE;

! Test comment

Q1:     QUADRUPOLE, L:=QL ;         ! check forward reference

Q1.F:   Q1,         K1=0.5 ;
Q1.D:   Q1,         K1:=-Q1.F->K1;  ! check attribute access
        Q1.F,       TILT=0.001;     ! check element update
SOL1:   SOLENOID,   L=0.5,
                    K1S=3.0;        ! check continuation
MULT1:  MULTIPOLE,  KNL={1.0, 1.0, 2.0, 6.0};
HK1:    HKICKER,    L=0,    KICK=0.001;
VK1:    VKICKER,    L=0,    KICK=-0.001;
HVK2:   KICKER      L=0     VKICK:=VKICK;  ! Check spaces as separators
BPM1:   MONITOR,    L=0.1;
BPM2:   VMONITOR;
RF:     RFCAVITY,   VOLT=5, FREQ=352.2, HARMON=31;

CELL.1: SEQUENCE,   L=4.0;
QFCELL: Q1.F,       AT=0.75;
        Q1.D,       AT=2.0, FROM=QFCELL;
ENDSEQUENCE;

RING.1: SEQUENCE,   L=LRING;
        CELL.1,     AT=2.0;
        SOL1,       AT=4.5;
        MULT1,      AT=4.8;
        HK1,        AT=4.8;
        VK1,        AT=4.8;
        HVK2,       AT=4.8;
        BPM1,       AT=4.9;
        BPM2,       AT=5;
        CELL.1,     AT=7.0;
ENDSEQUENCE;

VALUE, EMASS;

QL = sqrt(4)/2;          ! check arithmetic
const real lring = 9.0;
VKICK = 0.003;

TEST: TWISS, D=BETA;
"""

mad8_data = """
comment
test block comment
delimiter:      ";"
linecomment:    "!", "//"
endcomment

BEAM,   PARTICLE="positron", ENERGY=2.0;

! Test comment

Q1:     QUADRUPOLE, L:=QL ;         ! check forward reference

Q1.F:   Q1,         K1=0.5 ;
Q1.D:   Q1,         K1:=-Q1.F[K1];  ! check attribute access
        Q1.F,       TILT=0.001;     ! check element update
SOL1:   SOLENOID,   L=0.5,
                    KS=3.0;         ! check continuation
MULT1:  MULTIPOLE,  KNL={1.0, 1.0, 2.0, 6.0};
HK1:    HKICKER,    L=0,    KICK=0.001;
VK1:    VKICKER,    L=0,    KICK=-0.001;
HVK2:   KICKER,     L=0,    VKICK=VKICK;
BPM1:   MONITOR,    L=0.1;
BPM2:   VMONITOR;
RF:     RFCAVITY,   VOLT=5, FREQ=352.2, HARMON=31;

CELL.1: SEQUENCE,   L=4.0;
QFCELL: Q1.F,       AT=0.75;
        Q1.D,       AT=2.0, FROM=QFCELL;
ENDSEQUENCE;

RING.1: SEQUENCE,   L=LRING;
        CELL.1,     AT=2.0;
        SOL1,       AT=4.5;
        MULT1,      AT=4.8;
        HK1,        AT=4.8;
        VK1,        AT=4.8;
        HVK2,       AT=4.8;
        BPM1,       AT=4.9;
        BPM2,       AT=5;
        CELL.1,     AT=7.0;
ENDSEQUENCE;

VALUE, EMASS;

QL = sqrt(4)/2;          ! check arithmetic
LRING = 9.0;
VKICK = 0.003;

TEST: TWISS, D=BETA;
"""


elegant_data = """
! Test comment

DR1:    DRIFT,      L=0.25
DR2:    DRIFT,      L=1.0
DR3:    DRIFT,      L=0.75
DR4:    DRIFT,      L=0.05
Q1.F:   QUADRUPOLE, L="2 2 /", K1=0.5, TILT=0.001
Q1.D:   QUADRUPOLE, L=1.0, K1=-0.5 ! check attribute access
SOL1:   SOLENOID,   L=0.5, &
                    K1S=3.0        ! check continuation
MULT1:  MULTIPOLE,  KNL=6.0, ORDER=3
"HK:1": HKICKER,    L=0,    KICK=0.001
"VK:1": VKICKER,    L=0,    KICK=-0.001
HVK2:   KICKER,     L=0,    VKICK=0.003
BPM1:   MONITOR,    L=0.1
BPM2:   VMONITOR
RF:     RFCAVITY,   VOLT=5.0E6, FREQ=352.2E6, HARMON=31

CELL.1: LINE=(DR1, Q1.F, DR2, Q1.D, DR3)

RING.1: LINE=(CELL.1, DR1, SOL1, DR4, MULT1, &
              HK:1, VK:1, HVK2, DR4, BPM1, DR4, BPM2, CELL.1)
"""
# Global definitions for BaseParser

true = True
false = False


def _command1_parser(parser, argcount, argstr):
    return parser._argparser(
        argcount, argstr, bool_attr={"flag1", "FLAG2"}, str_attr={"title"}
    )


@set_argparser(_command1_parser)
def command1(**kwargs):
    """Sample command testing that arguments are as expected"""
    assert kwargs["title"] == "test parser"
    assert kwargs["flag1"] is True
    assert kwargs["FLAG2"] is False
    assert kwargs["arg1"] == 42.0
    assert kwargs["ARG2"] == 12.0
    assert kwargs["ARG3"] is True
    return "done"


# End of global definitions


@pytest.mark.parametrize(
    "cont, delimit, linecomm, data",
    [[None, ";", ("!", "//"), test_data1], ["\\", ";", "#", test_data2]],
)
def test_unordered_parser(cont, delimit, linecomm, data: str):
    class _Parser(UnorderedParser):
        _blockcomment = ("/*", "*/")
        _linecomment = linecomm
        _continuation = cont
        _delimiter = delimit

    parser = _Parser(globals())
    parser.parse_lines(data.splitlines())
    assert parser["label"] == "done"  # Make sure that "command1" was executed


@pytest.mark.parametrize(
    "Parser, data, energy, particle, radiate",
    [
        [MadxParser, madx_data, 6.0e9, "electron", True],
        [Mad8Parser, mad8_data, 2.0e9, "positron", False],
    ],
)
def test_madx_parser(Parser, data: str, energy: float, particle: str, radiate: bool):
    # fmt: off
    expected_ring_pos = np.array([
        0.0,    0.25,   1.25,   2.25,   3.25,   4.25,   4.75,   4.8,    4.8,
        4.8,    4.8,    4.8,    4.85,   4.9,    4.9,    4.95,   5.0,    5.0,
        5.25,   6.25,   7.25,   8.25,   9.0,
    ])
    expected_cell_pos = np.array([0.,   0.25,   1.25,   2.25,   3.25,   4.])
    # fmt: on
    parser = Parser()
    parser.parse_lines(data.splitlines())

    sequences = parser.sequences

    assert sequences == ["cell_1", "ring_1"]

    ring = parser.lattice(use="ring.1")
    assert len(ring) == 22
    assert ring.energy == energy
    assert str(ring.particle) == particle
    assert ring.circumference == 9.0
    assert ring.is_6d is radiate
    assert_close(ring.get_s_pos(), expected_ring_pos)

    cell = parser.lattice(use="cell.1")
    assert len(cell) == 5
    assert cell.energy == 2.0e09
    assert str(cell.particle) == "positron"
    assert cell.circumference == 4.0
    assert cell.is_6d is False
    assert_close(cell.get_s_pos(), expected_cell_pos)


def test_elegant_parser():
    # fmt: off
    expected_ring_pos = np.array([
        0.0, 0.25, 1.25, 2.25, 3.25, 4.0, 4.25, 4.75, 4.8, 4.8,
        4.8,  4.8,  4.8, 4.85,  4.9, 4.9, 4.95,  5.0, 5.0,
        5.25, 6.25, 7.25, 8.25, 9.0,
    ])
    expected_cell_pos = np.array([0., 0.25, 1.25, 2.25, 3.25, 4.])
    # fmt: on
    parser = ElegantParser()
    parser.parse_lines(elegant_data.splitlines())

    sequences = parser.sequences

    assert sequences == ["CELL_1", "RING_1"]

    ring = parser.lattice(use="ring.1", particle="electron", energy=6.0e9)
    assert len(ring) == 23
    assert ring.energy == 6.0e9
    assert str(ring.particle) == "electron"
    assert ring.circumference == 9.0
    assert ring.is_6d is False
    assert_close(ring.get_s_pos(), expected_ring_pos)

    cell = parser.lattice(use="cell.1", particle="positron", energy=2.0e9)
    assert len(cell) == 5
    assert cell.energy == 2.0e09
    assert str(cell.particle) == "positron"
    assert cell.circumference == 4.0
    assert cell.is_6d is False
    assert_close(cell.get_s_pos(), expected_cell_pos)
