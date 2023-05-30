from pdf_tonality import pdf_tonality as pdft
import numpy as np

def test_nll_tonality():

    # example 1:
    # an artificial pdf that contains only values of 1
    pitch_vec = np.repeat(12.0, 5)
    pitch_array = pdft.MIDI_to_transpose_array(pitch_vec)
    PDF = np.repeat(1, 1200)
    nll_array = pdft.nll_tonality(PDF, pitch_array)

    # log(1) = 0, so we expect an array of all zeros
    assert np.all(nll_array == 0)

    # example 2:
    # artificial "melody" with zeros for all notes, so all transpositions are
    # positive and first row is zeros
    pitch_vec = np.repeat(0.0, 5)
    pitch_array = pdft.MIDI_to_transpose_array(pitch_vec)
    PDF = np.arange(1, 1201)
    nll_array = pdft.nll_tonality(PDF, pitch_array)

    # log(1) = 0, so we expect an array of zeros in first row only, and only
    # non-zero after
    assert np.all(nll_array[0] == 0)
    assert np.all(nll_array[1:] != 0)

    # due to the way the transposition works in this case (always increasing),
    # and the way the artificial PDF was created (always increasing as well),
    # we expect only negative differences from one row to the next
    assert np.all(np.diff(nll_array, axis=0) < 0)