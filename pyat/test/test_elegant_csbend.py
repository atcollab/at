import numpy
import pytest

from at import element_pass, elements


@pytest.mark.parametrize(
    ("integration_order", "initial", "expected"),
    [
        (
            2,
            [8.0e-4, -3.0e-4, -6.0e-4, 4.0e-4, 7.0e-3, 0.0],
            [
                3.32012301394703701e-4,
                -4.14831730921788101e-4,
                -3.60191618748299119e-4,
                3.33108624500673542e-5,
                7.0e-3,
                5.64061899488829590e-4,
            ],
        ),
        (
            4,
            [-7.0e-4, 5.0e-4, 8.0e-4, -2.0e-4, -5.0e-3, 0.0],
            [
                -1.74379914127275595e-3,
                -1.89463221515014579e-3,
                2.70242759767550704e-4,
                -4.63045598364227886e-4,
                -5.0e-3,
                -9.50781118422522266e-4,
            ],
        ),
        (
            6,
            [8.0e-4, -3.0e-4, -6.0e-4, 4.0e-4, 7.0e-3, 0.0],
            [
                1.28397507302638911e-3,
                1.03676033774252781e-3,
                -3.64943348186758674e-4,
                2.67929594078719584e-5,
                7.0e-3,
                2.81346153117522106e-4,
            ],
        ),
    ],
)
def test_elegant_csbend_against_elegant(
    integration_order, initial, expected
):
    if integration_order == 4:
        length = 1.5
        angle = 0.7853981633974483
        k1 = -0.3
        k2 = 4.5
        slices = 4
    else:
        length = 1.2
        angle = 0.31
        k1 = 0.7
        k2 = -3.2
        slices = 5

    bend = elements.Dipole(
        "B",
        length,
        angle,
        PolynomB=[0.0, k1, k2 / 2.0],
        PassMethod="ElegantCsbendPass",
        NumIntSteps=slices,
        IntegrationOrder=integration_order,
        ExpansionOrder=4,
    )
    result = element_pass(bend, numpy.array(initial))

    numpy.testing.assert_allclose(result, expected, rtol=0.0, atol=2.0e-15)


def test_elegant_csbend_high_order_multipoles_against_elegant():
    initial = numpy.array(
        [1.1e-2, -2.0e-4, -7.0e-3, 3.0e-4, 3.0e-3, 0.0]
    )
    expected = numpy.array(
        [
            1.05839966826228003e-2,
            -1.07267156215775381e-3,
            -6.93054938968628957e-3,
            -8.45585350690351530e-5,
            2.99999999999989164e-3,
            1.40926800688029709e-3,
        ]
    )
    bend = elements.Dipole(
        "B",
        0.65,
        0.13,
        PolynomB=[
            0.0,
            0.12,
            -0.35,
            80.0,
            -8.0e3,
            6.0e5,
            -4.0e7,
            3.0e9,
            -2.0e11,
        ],
        PolynomA=[
            0.0,
            -0.03,
            0.08,
            -20.0,
            2.0e3,
            -1.5e5,
            1.0e7,
            -7.0e8,
            5.0e10,
        ],
        PassMethod="ElegantCsbendPass",
        NumIntSteps=8,
        IntegrationOrder=6,
        ExpansionOrder=0,
    )
    result = element_pass(bend, initial)

    numpy.testing.assert_allclose(
        result[:5], expected[:5], rtol=0.0, atol=2.0e-15
    )
    numpy.testing.assert_allclose(
        result[5], expected[5], rtol=0.0, atol=1.0e-9
    )


def test_elegant_csbend_auto_expansion_order():
    kwargs = dict(
        family_name="B",
        length=0.65,
        bending_angle=0.13,
        PolynomB=[0.0] * 8 + [-2.0e11],
        PassMethod="ElegantCsbendPass",
        NumIntSteps=8,
        IntegrationOrder=6,
    )
    initial = numpy.array(
        [1.1e-2, -2.0e-4, -7.0e-3, 3.0e-4, 3.0e-3, 0.0]
    )
    automatic = elements.Dipole(ExpansionOrder=0, **kwargs)
    explicit = elements.Dipole(ExpansionOrder=10, **kwargs)

    numpy.testing.assert_array_equal(
        element_pass(automatic, initial), element_pass(explicit, initial)
    )
