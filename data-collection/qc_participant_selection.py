#!/bin/python
import re
import pandas as pd

# Select the collection
collection_name = "GU_1"
# collection_name = "KKI_1"

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
# Rename columns
asd 	= asd.rename(columns={"participant_id": "bids_name"})
control = control.rename(columns={"participant_id": "bids_name"})

# Extract the target columns from MRIQC tsv files and save into a new dataframe
target_anat = group_anat[columns_anat]
target_func = group_func[columns_func]
# print(target_anat.head())
# print(target_func.head())

# Filter the MRIQC target data based on pre-defined quality control thresholds
target_anat = target_anat[(target_anat["cjv"] < cjv_max) &
						  (target_anat["cnr"] > cnr_min) &
						  (target_anat["efc"] < efc_max)].reset_index(drop=True)

target_func = target_func[(target_func["dvars_std"] < dvars_std_max) &
						  (target_func["fd_mean"] < fd_mean_max)].reset_index(drop=True)

# Find subjects that have sufficient anatomical and functional scan quality
## The bids_name column should only contain the participant ID
sub_label = r'sub-\d{5}'   # Regular expression for participant ID
for idx, row in target_anat.iterrows():
	sub_ID_anat  = re.findall(sub_label, row['bids_name'])
	target_anat.at[idx, 'bids_name'] = sub_ID_anat[0]

for idx, row in target_func.iterrows():
	sub_ID_func  = re.findall(sub_label, row['bids_name'])
	target_func.at[idx, 'bids_name'] = sub_ID_func[0]

common_rows = pd.merge(target_anat, target_func, on='bids_name')
print(f"Participants with sufficient quality scans: {len(common_rows)}")
common_asd 	   = pd.merge(common_rows, asd, 	on='bids_name')
common_control = pd.merge(common_rows, control, on='bids_name')
print(f"ASD participants: {len(common_asd)}")
print(f"TD participants: {len(common_control)}\n")

# Get the numbers of scan quality per diagnostic group
target_anat_asd = pd.merge(target_anat, asd, on='bids_name')
target_func_asd = pd.merge(target_func, asd, on='bids_name')
target_anat_control = pd.merge(target_anat, control, on='bids_name')
target_func_control = pd.merge(target_func, control, on='bids_name')
print(f"Anatomical scans with sufficient quality across subjects: {len(target_anat)}")
print(f"ASD participants: {len(target_anat_asd)}")
print(f"TD participants: {len(target_anat_control)}\n")
print(f"Functional scans with sufficient quality across subjects: {len(target_func)}")
print(f"ASD participants: {len(target_func_asd)}")
print(f"TD participants: {len(target_func_control)}\n")

# Save the count & subject list
with open(listdir + f"subList_{collection_name}_complete_quality_scan.txt", 'w') as file:
    file.write(f"Anatomical scans with sufficient quality: {len(target_anat)}\n")
    file.write(f"Functional scans with sufficient quality: {len(target_func)}\n")
    file.write("List of participant labels with sufficient scan quality in {}:\n".format(collection_name))
    for label in common_rows['bids_name']:
        file.write(label + '\n')

# Subset the dataframe with quality scans for subject id in the sublist: complete scan & questionnaire data
qc_subs    = common_rows[common_rows['bids_name'].isin(sublist['participant_id'])]
qc_asd 	   = qc_subs[qc_subs["bids_name"].isin(asd["bids_name"])]
qc_control = qc_subs[qc_subs["bids_name"].isin(control["bids_name"])]

print(f"Participants with complete quality scan data and complete questionnaire data: {len(qc_subs)}")
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