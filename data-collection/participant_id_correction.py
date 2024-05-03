#!/bin/python
import pandas as pd

# Set data directory
# datadir = "C:/LOCAL/PATH/GU_1/"   # Dataset from Georgetown University. Collection id: GU_1.
# datadir = "C:/LOCAL/PATH/KKI_1/"  # Dataset from Kennedy Krieger Institute. Collection id: KKI_1.
# datadir = "C:/LOCAL/PATH/OHSU_1/" # Dataset from Oregon Health and Science University. Collection id: OHSU_1.

# Load the spreadsheet that causes BIDS validation Error [Code 212] PARTICIPANT_ID_PATTERN
df = pd.read_csv(datadir + "participants.tsv", sep='\t', encoding='windows-1252')

# The required BIDS format of the values in the "participant_id" column: "sub-<subject_id>".
## Add "sub-" to each row in the "participant_id" column.
df['participant_id'] = 'sub-' + df['participant_id'].astype(str)

# Save the corrected spreadsheet
df.to_csv(datadir + "participants.tsv", sep='\t', index=False)
print("Resolved ERROR [CODE 212] PARTICIPANT_ID_PATTERN.")