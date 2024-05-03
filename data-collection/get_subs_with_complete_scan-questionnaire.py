#!/bin/python
import pandas as pd

## Load the participants.tsv files from the ABIDE II collections: subjects with complete scan data
collection_name = "GU_1"
# collection_name = "KKI_1"
# collection_name = "OHSU_1"
datadir 	  = "/LOCAL/PATH/"
participants  = pd.read_csv(datadir + f"{collection_name}/participants.tsv", sep='\t')

## Create dataframe with autistic participants and controls
### Diagnostic category is given in the dx_group column: 1=Autism, 2=Control
asd     = participants[participants['dx_group'] == 1].reset_index(drop=True)
control = participants[participants['dx_group'] == 2].reset_index(drop=True)

print(f"There are {len(asd)} ASD participants and {len(control)} controls in {collection_name} with complete scan data. Total: {len(participants)}.")

## Find autistic participants will complete target questionnaire data: ADI-R, SRS, and ADOS total scores
### Column names for ADI-R: adi_r_social_total_a, adi_r_verbal_total_bv, adi_r_nonverbal_total_bv, adi_r_rrb_total_c, adi_r_onset_total_d
### Column names for SRS: srs_total_raw, srs_total_t
### Column names for ADOS: ados_g_total, ados_2_total, ados_2_severity_total

### Drow rows (axis 0) containing missing values for a selected subset of column names. We don't exclude missing values for ADOS columns
### because a missing value does not necessarily mean that ADOS was not administered to the subject. It could mean that a different 
### ADOS version (ADOS G vs ADOS 2) was administered.

asd_complete = asd.dropna(axis=0, subset=['adi_r_social_total_a', 'adi_r_verbal_total_bv', 'adi_r_nonverbal_total_bv', 
						  'adi_r_rrb_total_c', 'adi_r_onset_total_d', 'srs_total_raw', 'srs_total_t']).reset_index(drop=True)
print(f"ASD participants with complete questionnaire data in {collection_name}: {len(asd_complete)}.")

## Stack the ASD participant datasets with complete questionnaire data onto the data of the controls.
complete_asd_hc_pheno = pd.concat([asd_complete, control], ignore_index=True)

## Save the data
complete_asd_hc_pheno.to_csv(f'ABIDEII-{collection_name}_complete_all_noqc.tsv', sep='\t', index=False)

## Create a text file with labels of the selected subjects: ASD participants with complete data (N = 33) and all controls (N = 56).
participant_list = complete_asd_hc_pheno['participant_id']
# print(participant_list)
with open(datadir + f'{collection_name}/subList_{collection_name}_complete_scan_questionnaire.txt', 'w') as f:
    for participant_id in participant_list:
        f.write(str(participant_id) + '\n')