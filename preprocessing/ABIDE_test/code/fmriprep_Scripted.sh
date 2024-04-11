#!/bin/bash
#Template provided by Daniel Levitas of Indiana University
#Edits by Andrew Jahn, University of Michigan, 07.22.2020
#Edits by Monika Utrosa Skerjanec, ENS-PSL, 31.03.2024
subject=$1

#User inputs:
bids_root_dir=$HOME/Documents/ABIDE_test/                     #Where is the data?
subj=$subject
echo $subj                          			                    #Which subject to analyze?
nthreads=4                          			                    #How many processors to use?
mem=20 #gb                                                    #How much memory to use in GB?
container=docker #docker or singularity

#Begin:

#Convert virtual memory from gb to mb
mem=`echo "${mem//[!0-9]/}"`       #remove gb at end, so fMRIprep can read without errors
mem_mb=`echo $(((mem*1000)-5000))` #reduce some memory for buffer space during pre-processing

#export TEMPLATEFLOW_HOME=$HOME/.cache/templateflow
#Where is the location of the Freesurfer license?
export FS_LICENSE=$HOME/Documents/ABIDE_test/derivatives/license.txt

#Run fmriprep
fmriprep-docker $bids_root_dir $bids_root_dir/derivatives participant \
  --participant-label $subj \
  --skip-bids-validation \
  --md-only-boilerplate \
  --fs-license-file $HOME/Documents/ABIDE_test/derivatives/license.txt \
  --fs-no-reconall \
  --skull-strip-t1w auto \
  --output-spaces MNI152NLin6Asym:res-2 \
  --nthreads $nthreads \
  --stop-on-first-crash \
  --mem_mb $mem_mb \
  -w $HOME 
  echo "diocane2"