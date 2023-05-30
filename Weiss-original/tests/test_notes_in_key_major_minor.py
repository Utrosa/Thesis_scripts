from pdf_tonality import pdf_tonality as pdft
import numpy as np
import pytest

def test_notes_in_key_major_minor():


    # example 1:
    # here we expect that 8 of 13 notes will be in key, and the key will be 'C'
    # in uppercase
    some_offset = 0.15
    pitch_vec = np.arange(0,13) + some_offset
    tonic = 0
    tonic_offset = some_offset
    mode = 'major'

    prop_in_key, key = pdft.notes_in_key_major_minor(
        pitch_vec = pitch_vec,
        tonic = tonic,
        tonic_offset = tonic_offset,
        mode = mode,
    )

    # expectation is for the proportion of notes in key to be approximately
    # 8/13
    expected = pytest.approx(8/13, abs=0.001)
    actual = prop_in_key
    assert actual == expected

    # specifying a tonic of 0 and a mode of major should return 'C' in uppercase
    assert key == 'C'

    # example 2:
    # here we expect that 7 of 12 notes will be in key, and the key will be 'f'
    # in lowercase
    some_offset = 0.15
    pitch_vec = np.arange(0,12) + some_offset
    tonic = 5
    tonic_offset = some_offset
    mode = 'minor'

    prop_in_key, key = pdft.notes_in_key_major_minor(
        pitch_vec = pitch_vec,
        tonic = tonic,
        tonic_offset = tonic_offset,
        mode = mode,
    )

    # expectation is for the proportion of notes in key to be approximately
    # 7/12
    expected = pytest.approx(7/12, abs=0.001)
    actual = prop_in_key
    assert actual == expected

    # specifying a tonic of 0 and a mode of minor should return 'f' in lowercase
    assert key == 'f'