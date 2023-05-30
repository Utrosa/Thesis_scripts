from pdf_tonality import pdf_tonality as pdft
import numpy as np
import scipy.stats as stats
import pytest

def test_generate_random_pitches():

    # example 1
    pitch_vec = np.tile(np.arange(10,30), 1000000)
    strategy = 'flat_octave'
    pitch_vec_rand = pdft.generate_random_pitches(pitch_vec, strategy)

    # expectation that the average pitch will be around 6 (half way from 0-12)
    actual = np.mean(pitch_vec_rand)
    expected = pytest.approx(6.0, abs=0.1)
    assert actual == expected

    # example 2
    pitch_vec = np.tile(
        np.concatenate([np.tile(10,10),np.tile(30,10)]),
        1000000
    )
    strategy = 'flat_range'
    pitch_vec_rand = pdft.generate_random_pitches(pitch_vec, strategy)

    # expectation that the average pitch will be around 20 (half way from 10-30)
    actual = np.mean(pitch_vec_rand)
    expected = pytest.approx(20.0, abs=0.1)
    assert actual == expected

    # expectation that the std pitch will be around np.sqrt(((max-min)**2)/12)
    actual = np.std(pitch_vec_rand)
    expected = pytest.approx(np.sqrt(((30-10)**2)/12), abs=0.1)
    assert actual == expected

    # example 3
    low_note = 1
    high_note = 100
    mu = 50
    sigma = 10
    size_notes = 1000000
    pitch_vec = np.round(
        stats.truncnorm.rvs(
            (low_note - mu) / sigma,
            (high_note - mu) / sigma,
            loc = mu,
            scale = sigma,
            size = size_notes
        ), decimals = 2
    )
    strategy = 'gaussian_range'
    pitch_vec_rand = pdft.generate_random_pitches(pitch_vec, strategy)

    # expectation that the average pitch will match the input
    actual = np.mean(pitch_vec_rand)
    expected = pytest.approx(np.mean(pitch_vec), abs=0.1)
    assert actual == expected

    # expectation that the std pitch will be match the input
    # note - this will only work when mu / sigma keep the curve far enough from
    # the boundaries (at least a few standard deviations)
    actual = np.std(pitch_vec_rand)
    expected = pytest.approx(np.std(pitch_vec), abs=0.1)
    assert actual == expected