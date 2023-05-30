from pdf_tonality import pdf_tonality as pdft
import numpy as np
from scipy import signal
import pytest

def test_PDF_null_notes():

    PDF, intervals, weights = pdft.build_diatonic_PDF('major')

    # expectation for all points in the PDF to sum to 1
    actual = np.sum(PDF)
    expected = pytest.approx(1.0, abs=0.000000001)
    assert actual == expected

    # expectation for peaks at certain values
    expected = np.array([200, 400, 500, 700, 900, 1100])
    actual = signal.find_peaks(PDF)[0]
    assert np.array_equal(actual, expected)

    # expectation for peaks at certain values
    PDF, _, _ = pdft.build_diatonic_PDF('minor')
    expected = np.array([200, 300, 500, 700, 800, 1000])
    actual = signal.find_peaks(PDF)[0]
    assert np.array_equal(actual, expected)