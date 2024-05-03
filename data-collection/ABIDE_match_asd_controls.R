#Remove all objects in the current workspace
rm(list=ls())

#Specify the collection name
# collection_name <- "GU_1"
# collection_name <- "KKI_1" # Exclude pilot subjects: 29273, 29274, 29275!
collection_name <- "OHSU_1"

#Specify the number of participants with complete quality scans and survey data
# subs_complete_data_all = 21 #GU
# subs_complete_data_all = 38 #KKI
subs_complete_data_all = 31 #OHSU: There are 32 subjects, but we need only 31*2

#### -------------------------------- START ------------------------------- ####
#Assign the file path & working directory
chemin= "C:/Users/monik/Documents/Paris/Academic Year (2023-24)/M2 Internship/Data-collection"
setwd(chemin)
input_dataset <- sprintf("ABIDEII-%s_complete_all_qc.tsv", collection_name)

#Load the packages
library(groupdata2)
library(tidyverse)
library(caret)
library(dplyr)

#### ------------------------------ GET DATA ------------------------------ ####
ABIDE=read.table(input_dataset,
                 sep='\t',
                 strip.white=TRUE,  # remove trailing/leading spaces from rows
                 header=TRUE, 
                 fill=TRUE)         # fill missing values with nan

#Replace NA values with "n/a" (BIDS requirement)
ABIDE[is.na(ABIDE)] <- "n/a"

#Subset the data based on subject type: autistic patients vs healthy controls
#For diagnostic group: 1 = Autism; 2 = Control.
subset_hs=(ABIDE$dx_group==2)
subset_asd=(ABIDE$dx_group==1)

#Create new data frames per subject type
ABIDE2_1=ABIDE[subset_hs,]
ABIDE2_3=ABIDE[subset_asd,]

#Stack the two new data sets on top of each other: controls at the top
ABIDE2=rbind(ABIDE2_1, ABIDE2_3)

#### ----------------------- PATIENT CONTROL SAMPLE ------------------------#### 
table(ABIDE2$sex, ABIDE2$dx_group) # Count the number of participants
set.seed(99)

#Sample ASD participants: 1=Autism; 2=Control
ASD=sample_n(ABIDE2[ABIDE2$dx_group==1,],subs_complete_data_all)
HC=sample_n(ABIDE2[ABIDE2$dx_group==2,],subs_complete_data_all)

test.data=rbind(ASD,HC)
table(test.data$dx_group)

#Age at the time of the functional scan
test.data$age_at_scan=as.numeric(test.data$age_at_scan)
subset_test_hs=test.data$dx_group==2
subset_test_asd=test.data$dx_group==1
t.test(test.data$age_at_scan [subset_test_hs],test.data$age_at_scan [subset_test_asd])

#Sex: 1=male; 2=female
table(test.data$sex, test.data$dx_group)
tablesex=table(test.data$sex,test.data$dx_group)
chisq.test(tablesex)

#Handedness: 1=right handed; 2=left handed; 3=mixed handed
table(test.data$handedness_category, test.data$dx_group)
tablemanuel=table(test.data$handedness_category,test.data$dx_group)
chisq.test(tablemanuel)

#Full Scale IQ (FIQ)
test.data$fiq=as.numeric(test.data$fiq)
t.test(test.data$fiq[subset_test_hs],test.data$fiq[subset_test_asd])

#Eyes status during resting scan: 1=open; 2=closed
unique_eye_status <- unique(test.data$eye_status_at_scan)
print(unique_eye_status)

#Are participants taking medication? 0=no; 1=yes
table(test.data$current_med_status, test.data$dx_group)
test.data$current_med_status=as.numeric(test.data$current_med_status)
t.test(test.data$current_med_status[subset_test_hs],test.data$current_med_status[subset_test_asd])

#What is the name of this medication?
print(test.data$current_medication_name[test.data$current_med_status == 1 & test.data$dx_group == 1])
print(test.data$current_medication_name[test.data$current_med_status == 1 & test.data$dx_group == 2])

#Off stimulants 24 hours prior to scan? 0=no; 1=yes
test.data$off_stimulants_at_scan[test.data$off_stimulants_at_scan == "n/a"] <- 0
test.data$off_stimulants_at_scan=as.numeric(test.data$off_stimulants_at_scan)
table(test.data$off_stimulants_at_scan, test.data$dx_group)
t.test(test.data$off_stimulants_at_scan[subset_test_hs],test.data$off_stimulants_at_scan[subset_test_asd])

#### ------------------ SAVE PARTICIPANTS WITH FULL DATA ------------------ ####
file_name <- sprintf("%s_matched_participants.tsv", collection_name)
write.table(file=file_name, test.data, row.names = FALSE, col.names = TRUE, 
            quote = FALSE, sep = '\t')
