from pdf_tonality import pdf_tonality as pdft
import numpy as np

def test_optimal_shift():

    # create a toy array with known characteristics
    nll_array = np.array(
        [
            [-10, -5, -1, -1],
            [-11, -6, -2, -2],
            [-12, -7, -3, -3], # minimal row if all note durations are equal
            [-11, -6, -2, -5], # minimal row if final note has greater duration
            [-10, -5, -1, -1],
        ]
    )

    # example 1
    # expectation for minimal index at row 2 for toy array, with equal pitch
    # durations (no input to pitch_durs)
    nll_means, nll_min_idx, nll_weighted = pdft.optimal_shift(nll_array)
    assert nll_min_idx == 2

    # expectation for nll_array == nll_weighted if no durations are provided
    assert np.array_equal(nll_array, nll_weighted)

    # expectation for means to equal the following array
    actual = nll_means
    expected = [-4.25, -5.25, -6.25, -6.0, -4.25]
    assert np.array_equal(actual, expected)

    # example 2
    # expectation for minimal index at row 3 for toy array, with final pitch
    # duration longer (hence more weighted) than the other three
    pitch_durs = np.array([1,1,1,5]) * .05
    _, nll_min_idx, _ = pdft.optimal_shift(nll_array, pitch_durs=pitch_durs)
    assert nll_min_idx == 3

    # example 3
    # the function should return all minimum values, if there are more than one.
    nll_array = np.array(
        [
            [-10, -5, -1, -1],
            [-11, -6, -2, -2],
            [-12, -7, -3, -3], # minimum row if all note durations are equal
            [-3, -3, -7, -12], # identical avg to min row if all durs equal
            [-10, -5, -1, -1],
        ]
    )

    # expectation for minimal indices at rows 2 and 3
    _, nll_min_idx, _ = pdft.optimal_shift(nll_array)
    actual = nll_min_idx
    expected = np.array([2,3])
    assert np.array_equal(actual, expected)