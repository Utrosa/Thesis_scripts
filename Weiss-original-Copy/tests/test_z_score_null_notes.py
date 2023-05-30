from pdf_tonality import pdf_tonality as pdft
import numpy as np
import pytest

def test_z_score_null_notes():

    # example 1:
    # a small null distribution of flat shape
    actual_PIKN = 0.55
    null_PIKN = np.arange(0.1, 1, 0.1)
    z_out = pdft.z_score_null_notes(actual_PIKN, null_PIKN, p_max=None)

    # for this example, we expect a percentile of 55.5
    assert z_out['pctl'] == pytest.approx(55.5, abs=.1)

    # for this example, we expect a z score of about 0.193 in the standard
    # method, and about 0.139 in the percentile method
    assert z_out['z'] == pytest.approx(0.193, abs=.001)
    assert z_out['z_pctl'] == pytest.approx(0.139, abs=.001)

    # example 2:
    # a larger null distribution of flat shape
    actual_PIKN = 0.55
    null_PIKN = np.arange(0.001, 1, 0.001)
    z_out = pdft.z_score_null_notes(actual_PIKN, null_PIKN, p_max=None)

    # for this example, we expect a percentile of 55.0
    assert z_out['pctl'] == pytest.approx(55.0, abs=.1)

    # for this example, we expect a z score of about 0.173 in the standard
    # method, and about 0.127 in the percentile method
    assert z_out['z'] == pytest.approx(0.173, abs=.001)
    assert z_out['z_pctl'] == pytest.approx(0.127, abs=.001)

    # example 3:
    # a larger null distribution of flat shape, score higher than p_max
    actual_PIKN = 1.0
    null_PIKN = np.arange(0.001, 1, 0.001)
    z_out = pdft.z_score_null_notes(actual_PIKN, null_PIKN, p_max=None)

    # for this example, we expect a percentile of 55.0
    assert z_out['pctl'] == pytest.approx(100.0, abs=.1)

    # for this example, we expect a z score of about 1.733 in the standard
    # method, and about 3.290 in the percentile method (the maximum limit from
    # p_max)
    assert z_out['z'] == pytest.approx(1.733, abs=.001)
    assert z_out['z_pctl'] == pytest.approx(3.290, abs=.001)

    # example 4:
    # setting a lower p_max reduces the z_pctl score in a predictable way
    actual_PIKN = 1.0
    null_PIKN = np.arange(0.001, 1, 0.001)
    p_max = .99 # smaller than the default based on null sequence length
    z_out = pdft.z_score_null_notes(actual_PIKN, null_PIKN, p_max=p_max)

    # for this example, we expect a z score of about 2.326 in the percentile
    # method (the maximum limit from p_max) but the standard method is
    # unchanged from the previous example
    assert z_out['z'] == pytest.approx(1.733, abs=.001)
    assert z_out['z_pctl'] == pytest.approx(2.326, abs=.001)