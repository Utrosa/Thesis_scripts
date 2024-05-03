#!/bin/python
import re
import pandas as pd

# Select the collection with multiple BOLD runs
collection_name = "OHSU_1"

#Set quality control thresholds
columns_anat    = ["bids_name", "cjv", "cnr", "efc"]    # Anatomical quality metrics
columns_func    = ["bids_name", "dvars_std", "fd_mean"] # Functional quality metrics
cjv_max 		= 0.8 # Coefficient of Joint Variation (lower values == better quality)
cnr_min 		= 1.2 # Contrast to Noise Ratio (higher values == better quality)
efc_max 		= 0.8 # Entropy Focus Criterion (lower values == better quality)
fd_mean_max 	= 0.5 # Framewise Displacement (lower values == better quality)
dvars_std_max 	= 1.5 # Normalized DVARS (lower values == better quality)

# Import the group MRIQC tsv files (participants with complete scan data, incomplete questionnaire data)
datadir = f"/LOCAL/PATH/{collection_name}/derivatives/mriqc/"
group_anat  = pd.read_csv(datadir + "group_T1w.tsv", sep='\t')
group_func  = pd.read_csv(datadir + "group_bold.tsv", sep='\t')

# Import the sublist with participants who have complete scan and complete questionnaire data
listdir = f"/LOCAL/PATH/{collection_name}/"
sublist = pd.read_csv(listdir + f"subList_{collection_name}_complete_scan_questionnaire.txt", header=None, names=['participant_id'])

# Import the phenotypic file to get stats for asd vs control
participants  = pd.read_csv(listdir + "participants.tsv", sep='\t')
asd    	      = participants[participants['dx_group'] == 1].reset_index(drop=True)
control       = participants[participants['dx_group'] == 2].reset_index(drop=True)
asd 		  = asd.rename(columns={"participant_id": "bids_name"})
control 	  = control.rename(columns={"participant_id": "bids_name"})


# Extract the target columns from MRIQC tsv files and save into a new dataframe
target_anat = group_anat[columns_anat]
target_func = group_func[columns_func]

# Filter the MRIQC target data based on pre-defined quality control thresholds
target_anat_qc = target_anat[(target_anat["cjv"] < cjv_max) &
						     (target_anat["cnr"] > cnr_min) &
						     (target_anat["efc"] < efc_max)].reset_index(drop=True)

target_func_qc = target_func[(target_func["dvars_std"] < dvars_std_max) &
						     (target_func["fd_mean"] < fd_mean_max)].reset_index(drop=True)

# How many quality BOLD runs remain per T1w run?
## Extract subject IDs from the quality anatomical data
sub_label   = r'sub-\d{5}'   # Regular expression for participant ID
run_label   = r'run-\d{1}'   # Regular expression for run number
sub_ID_anat = []
for idx, row in target_anat_qc.iterrows():
	sub_ID  = re.findall(sub_label, row['bids_name'])
	sub_ID_anat.extend(sub_ID)

scans_func = []
for idx, row in target_func_qc.iterrows():
	sub_ID  = re.findall(sub_label, row['bids_name'])
	run_ID  = re.findall(run_label, row['bids_name'])
	scan_ID = sub_ID[0] + '_' + run_ID[0]
	scans_func.append(scan_ID)

## Count how many times do the IDs repeat in the functional quality data. Save as dataframe.
occurrences = {sub_id: sum(sub_id in scan for scan in scans_func) for sub_id in sub_ID_anat}
count_df = pd.DataFrame(list(occurrences.items()), columns=['bids_name', 'BOLD_scans'])

print(f"Participants with sufficient quality scans: {len(count_df)}")
common_asd 	   = pd.merge(count_df, asd, 	 on='bids_name')
common_control = pd.merge(count_df, control, on='bids_name')
print(f"ASD participants: {len(common_asd)}")
print(f"TD participants: {len(common_control)}\n")

# Get the numbers of scan quality per diagnostic group
anat_asd     = pd.merge(count_df, asd,     on='bids_name')
anat_control = pd.merge(count_df, control, on='bids_name')
print(f"Anatomical scans with sufficient quality across subjects: {len(count_df)}")
print(f"ASD participants: {len(anat_asd)}")
print(f"TD  participants: {len(anat_control)}\n")

print(f"Functional scans with sufficient quality across subjects: {count_df['BOLD_scans'].sum()}")
bold_runs_asd  = count_df[count_df['bids_name'].isin(anat_asd['bids_name'])]
count_bold_asd = bold_runs_asd['BOLD_scans'].sum()
print(f"ASD participants: {count_bold_asd}")

bold_runs_control  = count_df[count_df['bids_name'].isin(anat_control['bids_name'])]
count_bold_control = bold_runs_control['BOLD_scans'].sum()
print(f"TD  participants: {count_bold_control}")

# Save a list of participants with quality scans 
with open(listdir + f"subList_{collection_name}_complete_quality_scan.txt", 'w') as file:
    file.write(f"Anatomical scans with sufficient quality: {len(count_df)}\n")
    file.write(f"Functional scans with sufficient quality: {count_df['BOLD_scans'].sum()}\n")
    file.write("List of participant labels with sufficient scan quality in {}:\n".format(collection_name))
    for label in count_df['bids_name']:
        file.write(label + '\n')

# Subset the dataframe with quality scans for subject id in the sublist: complete scan & questionnaire data
qc_subs    = count_df[count_df['bids_name'].isin(sublist['participant_id'])]
qc_asd 	   = qc_subs[qc_subs["bids_name"].isin(asd["bids_name"])]
qc_control = qc_subs[qc_subs["bids_name"].isin(control["bids_name"])]

print(f"\nParticipants with complete quality scan data and complete questionnaire data: {len(qc_subs)}")
print(f"ASD participants: {len(qc_asd)}")
print(f"TD participants: {len(qc_control)}")

# Save the subject list and a phenotypic file for matching
with open(listdir + f"subList_{collection_name}_complete_scan_questionnaire_quality.txt", 'w') as file:
    for label in qc_subs['bids_name']:
        file.write(label + '\n')

complete_qc_asd = asd[asd['bids_name'].isin(qc_asd['bids_name'])]
complete_qc_td  = control[control['bids_name'].isin(qc_control['bids_name'])]
complete_asd_hc_pheno = pd.concat([complete_qc_asd, complete_qc_td], ignore_index=True)
complete_asd_hc_pheno.to_csv(f'ABIDEII-{collection_name}_complete_all_qc.tsv', sep='\t', index=False)