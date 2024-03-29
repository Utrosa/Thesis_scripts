%%%  ERROR: Mask and Input volumes have different (x,y,z) size.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                        %
%        Parcellation of Preprocessed fMRI images using FSLmeants        %
%                                                                        %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Define paths -----------------------------------------------------------
addpath('/mnt/c/Users/monik/Documents/MATLAB/NIfTI_20140122/')
addpath('/home/moni/Documents/preprocessing/')
fprintf('\nBeginning the parcellation\n\n');
fsldir = '/usr/local/fsl/share/fsl/bin/';
setenv('FSLDIR','/usr/local/fsl');
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

dirhead        = pwd;
diratlas       = [dirhead '/atlas/'];
subjectdir     = [dirhead '/subject008/'];
atlas_filename = 'BN_Atlas_246_2mm.nii.gz';
filename       = 'f01_smoothed_WSL';
atlas          = 'fan246'; %Brainnetome atlas 2016 with 246 parcels

% Start the parcellation --------------------------------------------------
functional_parcellation_FSL(subjectdir, atlas_filename, filename, atlas)