
from pdf_tonality import pdf_tonality as pdft
import numpy as np

def test_PDF_null_notes():

    # happy birthday starting on MIDI 60 = C
    pitch_vec = np.array([60, 60, 62, 60, 65, 64, 60, 60, 62, 60, 67, 65, 60,
        60, 72, 69, 65, 64, 62, 70, 70, 69, 65, 67, 65], dtype='float')

    # add some random pitch error
    pitch_vec += np.random.normal(loc=0.0, scale=0.10, size=len(pitch_vec))

    # make one note just very slightly above zero, but less than 1 cent (0.01)
    # in MIDI notation, which would be rare but possible
    pitch_vec[0] = 0.001

    # calculate the transposition array
    pitch_array = pdft.MIDI_to_transpose_array(pitch_vec)

    # expectation for only differences of -1199 or 1 in the array from one row
    # to the next
    actual = np.unique(np.diff(pitch_array, axis=0))
    expected = np.array([-1199,1])
    assert np.array_equal(actual, expected)

    # expectation for the count of "-1199" to match the number of notes that
    # exceed or equal mod(n,12) of 0.01 (1 cent), since all notes of 1 cent or
    # more above middle C must wrap around exactly once
    actual = np.sum(np.diff(pitch_array, axis=0) == -1199)
    expected = np.sum(np.mod(pitch_vec, 12) >= 0.01)
    assert actual == expected