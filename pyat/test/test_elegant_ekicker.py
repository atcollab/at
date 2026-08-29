import numpy
import pytest

from at import element_pass, elements


@pytest.mark.parametrize(
    ("length", "kick_angle", "expected"),
    [
        (
            0.47,
            [2.3e-3, 0.0],
            [
                -6.83137897648110858e-3,
                -1.67002135572810739e-2,
                5.57809835306210134e-3,
                1.40005398036658936e-2,
                3.19999999999998064e-2,
                -1.88629997420554927e-3,
            ],
        ),
        (
            0.47,
            [0.0, -1.7e-3],
            [
                -7.35673935487284134e-3,
                -1.90003988184316921e-2,
                5.19042619321095043e-3,
                1.23002888711403781e-2,
                3.19999999999998064e-2,
                -1.88206340276794856e-3,
            ],
        ),
        (
            0.47,
            [2.3e-3, -1.7e-3],
            [
                -6.83124613115438078e-3,
                -1.67000065902711062e-2,
                5.19024607466582640e-3,
                1.23000036418330280e-2,
                3.19999999999998064e-2,
                -1.89135114494097751e-3,
            ],
        ),
        (
            0.0,
            [-3.1e-3, 2.2e-3],
            [
                1.29999999999999994e-3,
                -2.20990688769131458e-2,
                -7.99999999999999930e-4,
                1.61992887213782237e-2,
                3.19999999999998064e-2,
                -2.00000000000000004e-3,
            ],
        ),
    ],
)
def test_elegant_ekicker_against_elegant(length, kick_angle, expected):
    corrector = elements.Corrector(
        "C", length, kick_angle, PassMethod="ElegantEkickerPass"
    )
    initial = numpy.array(
        [1.3e-3, -1.9e-2, -8.0e-4, 1.4e-2, 3.2e-2, -2.0e-3]
    )

    result = element_pass(corrector, initial)

    numpy.testing.assert_allclose(result, expected, rtol=0.0, atol=3.0e-16)


def test_elegant_ekicker_zero_kick_matches_elegant_drift():
    corrector = elements.Corrector(
        "C", 0.63, [0.0, 0.0], PassMethod="ElegantEkickerPass"
    )
    drift = elements.Drift("D", 0.63, PassMethod="ElegantDriftPass")
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

    corrector_result = element_pass(corrector, initial.copy(order="F"))
    drift_result = element_pass(drift, initial.copy(order="F"))

    numpy.testing.assert_allclose(
        corrector_result, drift_result, rtol=0.0, atol=1.0e-17
    )
