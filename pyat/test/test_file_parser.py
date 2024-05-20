import pytest
from math import sqrt  # noqa: F401
from at import UnorderedParser

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
def test_unordered_parser1(delimiter, linecomment, data):
    parser = UnorderedParser(
        globals(),
        blockcomment=("/*", "*/"),
        linecomment=linecomment,
        delimiter=delimiter,
    )
    parser.parse_lines(test_data[data].splitlines())
    assert parser["label"] == "done"
