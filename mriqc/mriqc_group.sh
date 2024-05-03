#!/bin/bash

# Define inputs
nthreads=2
mem=10 #gb

# Run MRIQC
echo ""
echo "Running MRIQC on group level"
echo ""

# Replace LOCAL/PATH/abidedata with local path to your data folder
docker run -it --rm -v /LOCAL/PATH/abidedata:/data:ro \
	-v /LOCAL/PATH/abidedata/derivatives/mriqc:/out \
	nipreps/mriqc:24.0.0 /data /out \
	group \
	--n_proc $nthreads \
	--mem_gb $mem \
	--float32 \
	--no-sub
echo "DONE"