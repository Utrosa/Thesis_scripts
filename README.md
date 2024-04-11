# Dynamic Functional Connectivity (FC) in Autism Spectrum Disorder (ASD)
This repository contains code related to my master's thesis project. In short, the project entails power calculation, preprocessing resting-state fMRI data code, static FC analysis, and dynamic FC analysis.

## fMRIPrep preprocessing

A note on the arguments of fMRIPrep: we are analyzing data from children aged 8 -14. Since an unresolved bug prevents us from registering the data to a pediatric brain template, we registered the data to the MNI152NLin6Asym atlas. The good side of using this template is that for resolution 2mm (2x2x2 voxels), the shape of the atlas matches the shape of the atlas used for brain parcellation (91x106x91). This means that no resampling has to be performed. The Human brainnetome atlas only need to be re-oriented from Left-Right to Right-Left.

Run the fMRIPrep in bash, using fmriprep-docker, on all subjects. fmriprep-docker is a wrapper that translates fMRIPrep commands to docker commands so that no mounting has to be done manually. Since we have access to a high-performance computer, we can run the pipeline over all subjects in parallel. This speeds up the preprocessing.

There are two scripts related to this. Run the parallel_fMRIPrep script in the bash terminal, to start the preprocessing. If any errors occur, they will be saved in the log (see Docker logs or logs stored per subject folder in your output directory). Note that the pipeline will stop at the first crash. No worries though - you can re-run fMRIPrep with the same parameters, and the preprocessing will continue where it stopped. The -w parameter allows for this.

Sometimes fMRIPrep can hang ... re-run with low memory.

The boilerplate summarizing the preprocessing steps and parameters can be found at <OUTPUT_PATH>/logs/CITATION.md.

## Additional Preprocessing in MATLAB
These steps are largely based on the preprocessing steps from Abbas, Bassil et al. (2019) and Xu et al. (2023). The additional preprocessing has to be done because fMRIPrep skips certain steps, necessary for a Quasi-Periodic-Pattern (QPP) analysis.

These steps are trimming, spatial smoothing, filtering, signal regression (global, white matter, and cerebrospinal fluid signals), brain parcellation, and z-scoring. To run this step, the FMRIB Software Library (FSL) has to be installed. Therefore, this can only be performed on Linux systems or potentially through Docker.

Note that we will be performing a whole-brain analysis. Therefore, we will use the Human Brainnetome Atlas (Fan et al., 2016) to parcellate all voxels into 246 brain areas (parcels). Then, the mean time series for each brain parcel is extracted for each timepoint and saved into a MATLAB matrix, which serves as input for the QPP analysis.