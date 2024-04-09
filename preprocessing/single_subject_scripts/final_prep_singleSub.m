% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
%                 Smoothing, Filtering, Signal Regression,                %
%                    Brain Parcellation, and Z-scoring                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Note, this script will only work on a Linux system (or WSL for Windows).
clear; close all; clc
fprintf('\nBeginning the additional preprocessing pipeline\n\n');

%% SETUP ------------------------------------------------------------------
% Add the path to Nifti toolbox, which is necessary to read nii.gz files.
% Add the path to the main directory, where the Atlas and ABIDE data 
% folders are located.
addpath("/USER/MATLAB/NIfTI_20140122/", ...
        "/USER/preprocessing/")

%%% Set the FSL environment and output filetype to '.nii.gz'
fsldir = '/usr/local/fsl/share/fsl/bin/';
setenv('FSLDIR','/usr/local/fsl');
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

%%% Find fMRIPrep preprocessed data
dirhead         = pwd;
dirdata         = dir([dirhead '/ABIDE_test/derivatives']);
subjectdir      = [dirhead '/ABIDE_test/derivatives/sub-01/'];
filename        = 'sub-01_task-rest_run-1_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz';
% ABIDE_test is the folder with test data. The folder has to be in BIDS.

%%% Find the atlas
diratlas       = [dirhead '/Atlas/'];
atlas          = 'fan246'; %Brainnetome atlas 2016 with 246 parcels
% atlas_filename = [diratlas 'BN_Atlas_246_2mm.nii.gz'];
% reorient_atlas(diratlas, atlas_filename, atlas);
atlas_filename = [diratlas 'reorient_fan246.nii.gz'];
% To ensure a matching orientation between atlas and functional images run
% the reorient_atlas() function. If this has already been done, you can 
% skip this step and set the value of the atlas_filename variable directly 
% to the name of the reoriented atlas.

k              = strfind(subjectdir,'sub');  
sub            = [subjectdir(k:k+5),' - ']; 
% This is simply helping with the text that will be displayed as the
% function is running. 'sub' is a string with the name of the subject that
% will be printed with every update (see the Anzar Abbas' preprocessing 
% pipeline from 2016).

warning('off','all');
% There are a few warnings that MATLAB outputs depending on the type of
% computer this function that is being run on. As far as my understanding
% goes, they are irrelevant, hence the function is temporarily turning off
% warnings so that the updates being printed as the code runs look clean
% (see the Anzar Abbas' preprocessing pipeline from 2016).

%% SET THE PARAMETERS -----------------------------------------------------
%%% Trimming parameters
target_duration = 314000; % 04:26 is the shortest scan duration; here in ms
scan_duration   = 400000;
TR              = 2.340;
% OJO: TR for KKI_1 EPI scans is 2340ms whereas for GU_1 it's 2530ms. TR
% should be in seconds, not miliseconds!

%%% Smoothing parameters
smoothing_sigma = 6/sqrt(8*log(2)); %FWHM 6mm
% Smoothing sigma is needed for spatial smoothing with FSL. Keep in mind,
% authors usually report the full width at half maximum (FWHM), which can
% be easily transformed to sigma.
% FWHM = sigma*sqrt(8*ln(2)) = sigma*2.354

%%% Other parameters
fil            = [0.01, 0.08]; % Filtering bandwidth
gsr            = 1;            % 1 = yes and 0 = no
wmcsfr         = 1;            % 1 = yes and 0 = no
WM_probseg     = 'anat/sub-01_space-MNI152NLin6Asym_res-2_label-WM_probseg.nii.gz';
CSF_probseg    = 'anat/sub-01_space-MNI152NLin6Asym_res-2_label-CSF_probseg.nii.gz';

%% START ----------------------------------------------------------------
% Note, in MATLAB, you can call a function without keyword assignments!

% Trimming the functional scan to specified duration
trimming(subjectdir, ...
         fsldir, ...
         filename, ...
         target_duration, ...
         scan_duration);

% Spatial smoothing
spatial_smoothing(subjectdir, ...
                  fsldir, ...
                  smoothing_sigma);

filename = 'f0_smoothed';
% Temporal filtering and signal regression
filtering_signalreg(subjectdir, ...
                    filename, ...
                    TR, ...
                    fil, ...
                    gsr, ...
                    wmcsfr, ...
                     WM_probseg, ...
                    CSF_probseg);

filename = 'f0_smoothed_fil_gsr_wmcsfr';
% Brain parcellation and Z-scoring
functional_parcellation_zscore_AFNI(subjectdir, ...
                                    atlas_filename, ...
                                    filename, ...
                                    atlas);

%% End --------------------------------------------------------------------
time = whatsthetime();
fprintf([time,' ', sub,'Finished the final preprocessing steps.\n'])
warning('on','all');