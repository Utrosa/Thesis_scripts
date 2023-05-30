import os
import numpy as np
from rich import print
from pdf_tonality import pdf_tonality as pdft

# EXAMPLE 2: Demonstration of z-scores with a real melody, and with added
# "noise" to the pitch values.

# create a set of renditions of "happy birthday" that vary in accuracy, from
# perfectly in tune to large errors, or even complete randomness. If run
# repeatedly, we expect the z-scores to follow a descending order, though the
# results will vary from attempt to attempt (null z sd=0.045, normally
# distributed; most scores are within z=0.1 of mean, but at the extremes, could
# be slightly further).

# happy birthday in F-maj, 60=C, and we expect the proportion in key = 1.0
in_tune = np.array([60, 60, 62, 60, 65, 64, 60, 60, 62, 60, 67, 65, 60, 60, 72,
    69, 65, 64, 62, 70, 70, 69, 65, 67, 65])

# expect proportion in key = 1.0 but key of G (2.15 semitones above F)
in_tune_transposed = in_tune + 2.15 # +2 semitones and +15 cents to all pitches

# expect proportion in key = 0.96
one_bad_note = np.copy(in_tune)
one_bad_note[0] = 61

# expect proportion in key <= in_tune
small_error = in_tune + np.random.normal(
    loc = 0.0,
    scale = 0.30, # small error
    size = len(in_tune)
    ).round(2)

# expect proportion in key < small_error (but can vary, it is more or less
# random)
large_error = in_tune + np.random.normal(
    loc = 0.0,
    scale = 1.00, # large error
    size = len(in_tune)
    ).round(2)

# expect proportion in key to be very random, but with a z-score centered around
# 0 if run repeatedly.
random = np.round(60 + np.random.rand(len(in_tune)) * 12, decimals = 2)

pitch_vecs = {
    'in_tune': in_tune,
    'in_tune_transposed': in_tune_transposed,
    'one_bad_note': one_bad_note,
    'small_error': small_error,
    'large_error': large_error,
    'random': random,
}

# Load the pitch durations !
pitch_durs = np.array([])

# loop through the melodies and print results to console
for key in pitch_vecs:

    # print info
    pitch_vec = pitch_vecs[key]
    print(f"\nMelody: {key}")

    # create an "actual" model
    actual_model = pdft.PDF_major_minor(pitch_vec,pitch_durs)

    # create a null distribution
    null_summary, _ = pdft.PDF_null_notes(pitch_vec,pitch_durs)

    # calculate a z-score of the proportion of in-key notes
    z_out = pdft.z_score_null_notes(
        actual_PIKN=actual_model['prop_in_key'],
        null_PIKN=null_summary['prop_in_key'],
    )

    # print the results
    print(f"... the closest key is \'{actual_model['key']}\' (uppercase is major)")
    print(f"... the proportion of notes in key is: {z_out['actual_PIKN']}")
    print(f"... the percentile of that score is: {z_out['pctl']}")
    print(f"... the z-score (standard method) is {z_out['z'].round(2)}")
    print(f"... the z-score (percentile method) is {z_out['z_pctl'].round(2)}")