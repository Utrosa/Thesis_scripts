#!/bin/python
import re
import os
import shutil
import numpy as np
import pandas as pd

# Set the data directory & ABIDE II Collection Name ------------------------------------------------------------
## Dataset from Georgetown University. Collection id: GU_1.
#local_path = "/LOCAL/PATH"
datadir = f"{local_path}/GU_1/"
collection_name = "Georgetown University"
BIDS_report = "GU_1_missing_data_BIDS_error_report.txt"
# target_folder = f"{local_path}/GU_75"

# Resolve Warning 1: [Code 38] INCONSISTENT_SUBJECTS ------------------------------------------------------------
## Import the participant spreadsheet from source folder to extract all participant IDs
data = pd.read_csv(datadir + "participants.tsv", sep='\t', encoding='windows-1252')
participant_list = data['participant_id']
print(f'There are {len(participant_list)} participants in {collection_name} ABIDE II data collection.')

## Extract the participant labels from BIDS error & warning report. The labels in this report correspond to 
## participants who have missing functional data (BOLD scan).
collection_name == "Georgetown University"
## Define the regular expression for subject label
sub_label_missing_data =  r'sub-\d{5}'

## Open the report & read it
subjects = []
with open(datadir + BIDS_report, 'r') as report:
	text = report.read()
	subjs = re.findall(sub_label_missing_data, text)
	subjects.append(subjs)

## Extract the unique participant labels
missing_subjects = np.unique(subjects)
no_subjects = len(missing_subjects)
print(f'There are {no_subjects} subjects with missing data.')

# Select the participants for our study
no_selected_subjs = len(participant_list) - len(missing_subjects)
print(f'The number of selected subject from {collection_name} ABIDE II data collection is: {no_selected_subjs}')

## Extract only participant with full data
participant_list = participant_list.to_frame()
selected_subjs = participant_list[~participant_list['participant_id'].isin(missing_subjects)]
selected_subjs = selected_subjs.reset_index()

## Create a text file with labels of the selected subjects
with open(datadir + 'subjList.txt', 'w') as f:
	for index, row in selected_subjs.iterrows():
		if index == (len(selected_subjs) - 1):
			f.write(str(row['participant_id']))
		else:
			f.write(str(row['participant_id']) + '\n')

## Create a new folder with the selected subjects --------------------------------------------------------------

### Read subject list
with open(datadir + 'subjList.txt', 'r') as sf:
	folders_to_move = sf.read().splitlines()

### Move folders
# count = 0
# for folder_name in folders_to_move:
#     source_path = os.path.join(datadir, folder_name)
#     target_path = os.path.join(target_folder, folder_name)
#     if os.path.exists(source_path):
#         shutil.move(source_path, target_path)
#         print(f"Moved {folder_name} to {target_folder}")
#         count += 1
#     else:
#         print(f"Folder {folder_name} not found in {datadir}")

# print(f"Total {count} folders moved to {target_folder}")