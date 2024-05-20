import pytest
from math import sqrt  # noqa: F401
from at import UnorderedParser, MadxParser

__all__ = []

true = True
false = False

test_data1 = """
/* test block comment
delimiter:      ";"
linecomment:    "!", "//"
end of comments */

V1 = 42.0 ;                  // Test line comment
V2 = 2*(sqrt(4) + 2**2) ;    // Test expression

label : command1, flag1,     ! Test continuation
-flag2, title="test parser", ! Test disabled flag
arg1=v1, arg2=v2, arg3=v3 ;  ! Test postponed definition

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
-flag2, title="test parser",\ # Test disabled flag
arg1=v1, arg2=v2, arg3=v3     # Test postponed definition

V3 = True
"""

test_data = dict(data1=test_data1, data2=test_data2)

madx_data = """
BEAM, PARTICLE='electron', RADIATE, ENERGY=6.0, SEQUENCE=RING;
BEAM, PARTICLE='positron', ENERGY=2.0;

Q1: QUADRUPOLE, L:=QL ;  ! check forward reference

QF1: Q1, K1=0.5 ;
QD1: Q1, K1=-QF1->K1;    ! check attribute access
QF1, TILT=0.001;         ! check element update 
SOL1: SOLENOID, L=0.5, K1S=3.0;
MULT1: MULTIPOLE, KNL={1.0, 1.0, 2.0, 6.0};
HK1: HKICKER, L=0, KICK=0.001;
VK1: VKICKER, L=0, KICK=-0.001;
HVK2: KICKER, L=0, VKICK=VKICK;
BPM1: MONITOR, L=0.1;
BPM2: VMONITOR;
RFCAVITY, VOLT=5, FREQ=352.2, HARMON=31, AT=0.0;

RING: SEQUENCE, L=9.0;
SOL1, AT=4.5;
MULT1, AT=4.8;
HK1, AT=4.8;
VK1, AT=4.8;
HVK2, AT=4.8;
BPM1, AT=4.9;
BPM2, AT=5;
ENDSEQUENCE;

VALUE, EMASS;

QL = sqrt(4)/2;          ! check arithmetic
VKICK = 0.003;
"""


def command1(**kwargs):
    """Sample command testing that arguments are as expected"""
    print(kwargs)
    assert kwargs["title"] == "test parser"
    assert kwargs["flag1"] is True
    assert kwargs["flag2"] is False
    assert kwargs["arg1"] == 42.0
    assert kwargs["arg2"] == 12.0
    assert kwargs["arg3"] is True
    return "done"


@pytest.mark.parametrize(
    "delimiter, linecomment, data", [[";", ("!", "//"), "data1"], [None, "#", "data2"]]
)
def test_unordered_parser(delimiter, linecomment, data):
    parser = UnorderedParser(
        globals(),
        blockcomment=("/*", "*/"),
        linecomment=linecomment,
        delimiter=delimiter,
    )
    parser.parse_lines(test_data[data].splitlines())
    assert parser["label"] == "done"


def test_madx_parser():
    parser = MadxParser()
    parser.parse_lines(madx_data.splitlines())
    ring = parser.lattice(use="ring")
