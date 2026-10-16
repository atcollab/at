import numpy

from at import element_pass, elements


def test_elegant_drift_against_elegant():
    drift = elements.Drift("D", 1.7, PassMethod="ElegantDriftPass")
    initial = numpy.array([1.0e-3, 0.23, -2.0e-3, -0.17, 0.08, 4.0e-3])

    result = element_pass(drift, initial)

    # Generated with elegant EDRIFT and converted from slopes to AT momenta.
    expected = numpy.array(
        [
            3.76441224231789107e-1,
            2.3e-1,
            -2.79500035301757210e-1,
            -1.7e-1,
            8.0e-2,
            6.69414007405751549e-2,
        ]
    )
    numpy.testing.assert_allclose(result, expected, rtol=0.0, atol=2.0e-15)


def test_elegant_drift_matches_exact_hamiltonian():
    elegant_drift = elements.Drift("D", -0.9, PassMethod="ElegantDriftPass")
    exact_drift = elements.Drift("D", -0.9, PassMethod="ExactDriftPass")
    initial = numpy.array(
        [
            [1.0e-3, -2.0e-3],
            [0.31, -0.27],
            [-3.0e-3, 4.0e-3],
            [-0.22, 0.19],
            [0.12, -0.07],
            [5.0e-3, -6.0e-3],
        ],
        order="F",
    )

    elegant_result = element_pass(elegant_drift, initial.copy(order="F"))
    exact_result = element_pass(exact_drift, initial.copy(order="F"))

    numpy.testing.assert_equal(elegant_result, exact_result)
