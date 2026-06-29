import numpy
import pytest

from at import element_pass, elements


@pytest.mark.parametrize(
    ("integration_order", "expected"),
    [
        (
            2,
            [
                -4.38688681054565099e-3,
                -1.80049069307838677e-2,
                3.28319849111700403e-3,
                1.29830961085751206e-2,
                2.70000000000001350e-2,
                -2.91350733736464237e-3,
            ],
        ),
        (
            4,
            [
                -4.38689132876458540e-3,
                -1.80049693676556655e-2,
                3.28319021497942192e-3,
                1.29829089874033844e-2,
                2.70000000000001350e-2,
                -2.91350736328984923e-3,
            ],
        ),
        (
            6,
            [
                -4.38689132933226065e-3,
                -1.80049693655000669e-2,
                3.28319021525791703e-3,
                1.29829089852319113e-2,
                2.70000000000001350e-2,
                -2.91350736326889464e-3,
            ],
        ),
    ],
)
def test_elegant_ksext_against_elegant(integration_order, expected):
    sext = elements.Sextupole(
        "S",
        0.37,
        6.4,
        PassMethod="ElegantKsextPass",
        NumIntSteps=8,
        IntegrationOrder=integration_order,
    )
    initial = numpy.array(
        [2.1e-3, -1.8e-2, -1.4e-3, 1.3e-2, 2.7e-2, -3.0e-3]
    )

    result = element_pass(sext, initial)

    numpy.testing.assert_allclose(result, expected, rtol=0.0, atol=2.0e-16)


def test_elegant_ksext_zero_strength_matches_elegant_drift():
    sext = elements.Sextupole(
        "S",
        0.9,
        0.0,
        PassMethod="ElegantKsextPass",
        NumIntSteps=6,
        IntegrationOrder=6,
    )
    drift = elements.Drift("D", 0.9, PassMethod="ElegantDriftPass")
    initial = numpy.array(
        [
            [1.0e-3, -2.0e-3],
            [0.19, -0.15],
            [-3.0e-3, 4.0e-3],
            [-0.12, 0.10],
            [0.07, -0.03],
            [5.0e-3, -6.0e-3],
        ],
        order="F",
    )

    sext_result = element_pass(sext, initial.copy(order="F"))
    drift_result = element_pass(drift, initial.copy(order="F"))

    numpy.testing.assert_allclose(
        sext_result, drift_result, rtol=0.0, atol=4.0e-16
    )


def test_elegant_ksext_scaling_and_steering_against_elegant():
    sext = elements.Sextupole(
        "S",
        0.31,
        -7.6,
        PassMethod="ElegantKsextPass",
        NumIntSteps=9,
        IntegrationOrder=4,
        FieldScaling=1.06,
        KickAngle=numpy.array([1.3e-4, -1.9e-4]),
    )
    initial = numpy.array(
        [-1.5e-3, 1.4e-2, 1.2e-3, -1.1e-2, -2.2e-2, 2.5e-3]
    )
    expected = numpy.array(
        [
            2.95914420653304320e-3,
            1.41321156116792686e-2,
            -2.31672171864316006e-3,
            -1.11814653620171886e-2,
            -2.19999999999997975e-2,
            2.55201455169607209e-3,
        ]
    )

    result = element_pass(sext, initial)

    numpy.testing.assert_allclose(result, expected, rtol=0.0, atol=3.0e-16)
