#!/bin/bash

# Define paths & other inputs
# Replace LOCAL/PATH/abidedata with your local path
bids_root_dir=$HOME/LOCAL/PATH/abidedata
subject=$1
subj=$subject
nthreads=2
mem=10 #gb

# Run MRIQC
echo ""
echo "Running MRIQC on participant $subj"
echo ""

docker run -it --rm -v /LOCAL/PATH/abidedata:/data:ro \
	-v /LOCAL/PATH/abidedata/derivatives/mriqc:/out \
	nipreps/mriqc:24.0.0 /data /out \
	participant \
	--participant-label $subj \
 	--n_proc $nthreads \
	--mem_gb $mem \
	--float32 \
	--no-sub \
	--work-dir $bids_root_dir/derivatives/mriqc
echo "DONE"