from pdf_tonality import pdf_tonality as pdft
import numpy as np

def test_PDF_major_minor():

    # example 1: pitches expected to generate a *major* response with tonic = 0
    pitch_vec = np.array([0,0,2,4,5,7,9,11]) + 0.15
    pitch_durs = np.ones_like(pitch_vec)
    model_out = pdft.PDF_major_minor(pitch_vec, pitch_durs)

    assert model_out['params'] == 'major'
    assert model_out['tonic'] == 0
    assert model_out['tonic_offset'] == 0.15

    # example 2: pitches expected to generate a *minor* response with tonic = 0
    pitch_vec = np.array([0,0,2,3,5,7,8,10]) + 0.15
    pitch_durs = np.ones_like(pitch_vec)
    model_out = pdft.PDF_major_minor(pitch_vec, pitch_durs)

    assert model_out['params'] == 'minor'
    assert model_out['tonic'] == 0
    assert model_out['tonic_offset'] == 0.15