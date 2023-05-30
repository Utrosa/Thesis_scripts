from pdf_tonality import pdf_tonality as pdft
import numpy as np
import pytest

def test_run_PDF():

    # example 1:
    # expect two pitch vectors with a non-octave offset in semitones (int) to
    # result in the same nll and tonic offset, but different tonics
    pitch_vec1 = np.array([0,0,0,2,4,5,7,9,11,3]) + 0.3
    pitch_vec2 = pitch_vec1 + np.random.randint(1,11)
    model1 = pdft.run_PDF('major', pitch_vec1)
    model2 = pdft.run_PDF('major', pitch_vec2)

    print(model1)

    # expect the key to be C=0 in model1, and 9 of 10 notes to be in key
    assert model1['key'] == 'C'
    assert model1['prop_in_key'] == pytest.approx(9/10, abs=0.001)

    # expectation for best fit nll to match
    assert model1['nll'] == model2['nll']

    # expectation for tonic offset to match
    assert model1['tonic_offset'] == model2['tonic_offset']

    # expectation for tonics to differ
    assert model1['tonic'] != model2['tonic']