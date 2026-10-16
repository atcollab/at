import numpy
import pytest

from at import element_pass, elements


@pytest.mark.parametrize(
    ("integration_order", "expected"),
    [
        (
            2,
            [
                -1.31706973608012399e-2,
                -1.44342229811415378e-2,
                1.28155817114790904e-2,
                2.43117490896851661e-2,
                3.09999999999999165e-2,
                -3.72743651577958548e-3,
            ],
        ),
        (
            4,
            [
                -1.31955436838391300e-2,
                -1.44451727978595849e-2,
                1.27957979833768897e-2,
                2.43217496373718221e-2,
                3.09999999999999165e-2,
                -3.72753177981045283e-3,
            ],
        ),
        (
            6,
            [
                -1.31955956651026007e-2,
                -1.44447756892513903e-2,
                1.27960183622720693e-2,
                2.43221328016782687e-2,
                3.09999999999999165e-2,
                -3.72752267828957342e-3,
            ],
        ),
    ],
)
def test_elegant_kquad_against_elegant(integration_order, expected):
    quad = elements.Quadrupole(
        "Q",
        0.73,
        1.8,
        PassMethod="ElegantKquadPass",
        NumIntSteps=7,
        IntegrationOrder=integration_order,
    )
    initial = numpy.array(
        [1.2e-3, -2.3e-2, -0.8e-3, 1.7e-2, 3.1e-2, -4.0e-3]
    )

    result = element_pass(quad, initial)

    numpy.testing.assert_allclose(result, expected, rtol=0.0, atol=2.0e-16)


def test_elegant_kquad_zero_strength_matches_elegant_drift():
    quad = elements.Quadrupole(
        "Q",
        1.1,
        0.0,
        PassMethod="ElegantKquadPass",
        NumIntSteps=5,
        IntegrationOrder=6,
    )
    drift = elements.Drift("D", 1.1, PassMethod="ElegantDriftPass")
    initial = numpy.array(
        [
            [1.0e-3, -2.0e-3],
            [0.21, -0.17],
            [-3.0e-3, 4.0e-3],
            [-0.13, 0.11],
            [0.08, -0.04],
            [5.0e-3, -6.0e-3],
        ],
        order="F",
    )

    quad_result = element_pass(quad, initial.copy(order="F"))
    drift_result = element_pass(drift, initial.copy(order="F"))

    numpy.testing.assert_allclose(
        quad_result, drift_result, rtol=0.0, atol=4.0e-16
    )


def test_elegant_kquad_scaling_and_steering_against_elegant():
    quad = elements.Quadrupole(
        "Q",
        0.41,
        -2.3,
        PassMethod="ElegantKquadPass",
        NumIntSteps=9,
        IntegrationOrder=4,
        FieldScaling=1.07,
        KickAngle=numpy.array([1.4e-4, -2.1e-4]),
    )
    initial = numpy.array(
        [-7.0e-4, 1.1e-2, 9.0e-4, -1.5e-2, -2.4e-2, 3.0e-3]
    )
    expected = numpy.array(
        [
            4.13204920892837682e-3,
            1.28080552318562899e-2,
            -5.19310730967660703e-3,
            -1.29719477227512123e-2,
            -2.39999999999999103e-2,
            3.07389403864655013e-3,
        ]
    )

    result = element_pass(quad, initial)

    numpy.testing.assert_allclose(result, expected, rtol=0.0, atol=2.0e-16)


def test_elegant_kquad_fringe_against_elegant():
    fringe_int_m = numpy.array(
        [
            -9.7421568360565817e-3,
            2.6458577269830992e-4,
            -1.5177531216328479e-5,
            1.2817867420570844e-6,
            1.3456484106907799e-6,
        ]
    )
    fringe_int_p = numpy.array(
        [
            9.7421568360565817e-3,
            2.6458577269830992e-4,
            1.5177531216328479e-5,
            1.2817867420570844e-6,
            1.3456484106907799e-6,
        ]
    )
    quad = elements.Quadrupole(
        "Q",
        0.73,
        1.8,
        PassMethod="ElegantKquadPass",
        NumIntSteps=7,
        IntegrationOrder=6,
        FringeQuadEntrance=3,
        FringeQuadExit=3,
        fringeIntM0=fringe_int_m,
        fringeIntP0=fringe_int_p,
    )
    initial = numpy.array(
        [1.2e-3, -2.3e-2, -0.8e-3, 1.7e-2, 3.1e-2, -4.0e-3]
    )
    expected = numpy.array(
        [
            -1.3220742230563189e-2,
            -1.4440055732329996e-2,
            1.2770610486806150e-2,
            2.4317673528403555e-2,
            3.1e-2,
            -3.727512871431733e-3,
        ]
    )

    result = element_pass(quad, initial)

    numpy.testing.assert_allclose(result, expected, rtol=0.0, atol=2.0e-16)
