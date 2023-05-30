import os
import numpy as np
import pandas as pd
from pdf_tonality import pdf_tonality as pdft

# EXAMPLE 1: Demonstration of the PDF_major_minor() function to identify the key
# and proportion of in-key notes for a given set of pitches, followed by the
# PDF_null_notes() function, which creates a null distribution of proportion of
# in-key notes by repeatedly running PDF_major_minor() with random pitches but
# the same durations as the initial melody. To save the output files that get
# generated in the process, specify a path below and comment out the line that
# resets the path ("dir_out_detailed = None").

# if desired, set a destination for the full output to be written, and toggle
# write_output to "True"
write_output = True
dir_out = "C:/Users/monik/Documents/Paris/Internship/RESULTS/Weiss"
participant_ID = ''
condition = 'IMPROV'

# some setup
if write_output == True:
    # check that the path provided exists
    assert os.path.exists(dir_out), "Provide a valid path"
    # WARNING: if dir_out_detailed does not exist, it will be created!
    dir_out_detailed = os.path.join(dir_out, f'pdft_pilotdata_{participant_ID}_{condition}/')
    save_method = 'summary'
    plot_model = True
else:
    #dir_out_detailed = None
    save_method = 'summary'
    plot_model = False

# Create a set of pitches: load clean IMPROV data.
df = pd.read_csv(f"C:/Users/monik/Documents/Paris/Internship/Data-Collection/Weiss-Like/Weiss_data_{participant_ID}_{condition}.csv")

# Select the MIDI NOTE column without NaNs & Turn the MIDI notes into a numpy array.
MIDI_pilot = df['NOTE'].dropna().to_numpy()
DURATION_pilot = df['DURATION'].dropna().to_numpy()

pitch_vec = np.array(MIDI_pilot)
pitch_durs = np.array(DURATION_pilot)

# create an "actual" model, i.e., using the notes specified -----------------------------------------------------------
actual_model = pdft.PDF_major_minor(
    pitch_vec=pitch_vec,
    pitch_durs=pitch_durs,
    dir_out_detailed=dir_out_detailed,
    save_method=save_method,
    plot_model=plot_model,
)

# now, create a null distribution based on the actual durations, but generated
# from random pitch values on each of n iterations

# plotting each model in this case is incredibly slow and should be avoided

# saving the full output of each model slows the progression by ~5-10x and
# should be in a new directory, hence os.path.join(). The summary output is not
# too much slower per run.

# WARNING as with the code above, this output folder will be created!

dir_out_detailed_null = None
if dir_out_detailed is not None:
    dir_out_detailed_null = os.path.join(dir_out_detailed,'null_models')

n = 1000 # Number of iterations for random pitches ---------------------------------------------------------------------
null_summary, null_pitches = pdft.PDF_null_notes(
    pitch_vec=pitch_vec,
    pitch_durs=pitch_durs,
    strategy='gaussian_range',
    n=n,
    dir_out_detailed=dir_out_detailed_null,
    save_method='summary',
    plot_model=False,
)

# save the output if desired
if dir_out_detailed is not None:
    null_summary.to_csv(
        os.path.join(dir_out_detailed, f'null_models_summary_{participant_ID}_{condition}.csv'),
        index=False,
    )
    null_pitches.to_csv(
        os.path.join(dir_out_detailed, f'null_models_pitches_{participant_ID}_{condition}.csv'),
        index=False,
    )

# finally, we can calculate a z-score of the proportion of in-key notes and ------------------------------------------------------
# display the results
z_out = pdft.z_score_null_notes(
    actual_PIKN=actual_model['prop_in_key'],
    null_PIKN=null_summary['prop_in_key'],
)

print(f"\nThe proportion of notes for {participant_ID} in key is: {z_out['actual_PIKN']}, percentile {z_out['pctl']}")
print(f"\nThe z-score (standard method) for {participant_ID} is {z_out['z']}")
print(f"\nThe z-score (percentile method) for {participant_ID} is {z_out['z_pctl']}")

# show a plot of the z-score distribution
pdft.plot_PIKN_distribution(
    value=z_out['actual_PIKN'],
    distribution=z_out['null_PIKN']
)