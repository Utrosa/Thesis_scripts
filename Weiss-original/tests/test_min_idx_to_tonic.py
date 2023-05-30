from pdf_tonality import pdf_tonality as pdft
import numpy as np

def test_min_idx_to_tonic():

    # example 1: a tonic 10 cents below C=0 would have a minimum index of 10
    nll_min_idx = 10
    tonic, tonic_offset = pdft.min_idx_to_tonic(nll_min_idx)

    assert tonic == 0
    assert tonic_offset == -0.10

    # example 2: a tonic 50 cents above D=2 would have a minimum index of 950
    nll_min_idx = 950
    tonic, tonic_offset = pdft.min_idx_to_tonic(nll_min_idx)

    assert tonic == 2
    assert tonic_offset == 0.50

    # example 3: a tonic 49 cents below D=2 would have a minimum index of 1049
    nll_min_idx = 1049
    tonic, tonic_offset = pdft.min_idx_to_tonic(nll_min_idx)

    assert tonic == 2
    assert tonic_offset == -0.49