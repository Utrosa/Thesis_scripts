import os
import numpy as np
import pandas as pd
from rich import print
from pdf_tonality import pdf_tonality as pdft

# Load clean IMPROV data.
df = pd.read_csv("C:/Users/monik/Documents/Paris/Internship/Analyses/Pilot-Data/Big-Sets/Clean/spectral_data_Clean_ALL_IMPROV_pilots.csv")

# Select the MIDI NOTE column without NaNs & Turn the MIDI notes into a numpy array.
MIDI_pilot = df['MIDI_NOTE'].dropna().to_numpy()

# Specify the sung improvsiations to analyze
in_tune = np.array(MIDI_pilot)

# Get a random melody & Expect proportion in key to be very random, 
# but with a z-score centered around 0 if run repeatedly.
random = np.round(60 + np.random.rand(len(in_tune)) * 12, decimals = 2)
    
pitch_vecs = {
    'in_tune': in_tune,
    'random': random,
}

# Get DURATIONS of each MIDI note
pitch_durs = np.ones(len(MIDI_pilot)) ## Assume equal duration for all notes for now !!!!

# loop through the melodies and print results to console
for key in pitch_vecs:

    # print info
    pitch_vec = pitch_vecs[key]
    print(f"\nMelody: {key}")

    # create an "actual" model
    actual_model = pdft.PDF_major_minor(pitch_vec,pitch_durs)

    # create a null distribution
    null_summary, _ = pdft.PDF_null_notes(pitch_vec,pitch_durs, strategy='gaussian_range', n=1000,
    dir_out_detailed=None, save_method='summary', plot_model=False)

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