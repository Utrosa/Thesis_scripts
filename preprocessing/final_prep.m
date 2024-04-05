% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
%                 Smoothing, Filtering, Signal Regression,                %
%                    Brain Parcellation, and Z-scoring                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
clear; close all; clc
fprintf('\nBeginning the additional preprocessing pipeline\n\n');

%% SETUP ------------------------------------------------------------------
addpath("/mnt/c/Users/monik/Documents/MATLAB/NIfTI_20140122/", ...
        "/home/moni/Documents/preprocessing/")

% Set the FSL environment and output filetype to '.nii.gz'
fsldir = '/usr/local/fsl/share/fsl/bin/';
setenv('FSLDIR','/usr/local/fsl');
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

% Find fMRIPrep preprocessed data
dirhead         = pwd;
dirdata         = dir([dirhead '/ABIDE_test/derivatives']);
subjectdir      = [dirhead '/ABIDE_test/derivatives/sub-01/'];
filename        = 'sub-01_task-rest_run-1_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz';

% Find the atlas
diratlas       = [dirhead '/Atlas/'];
atlas          = 'fan246'; %Brainnetome atlas 2016 with 246 parcels
atlas_filename = [diratlas 'BN_Atlas_246_2mm.nii.gz'];

k              = strfind(subjectdir,'sub');  
sub            = [subjectdir(k:k+5),' - ']; 
% This is simply helping with the text that will be displayed as the
% function is running. 'sub' is a string with the name of the subject that
% will be printed with every update (see Anzar Abbas, 2016).

warning('off','all');
% There are a few warnings that MATLAB outputs depending on the type of
% computer this function that is being run on. As far as my understanding
% goes, they are irrelevant, hence the function is temporarily turning off
% warnings so that the updates being printed as the code runs look clean
% (see Anzar Abbas, 2016).

%% SET THE PARAMETERS -----------------------------------------------------

smoothing_sigma = 6/sqrt(8*log(2)); %FWHM 6mm
% Smoothing sigma is needed for spatial smoothing with FSL. Keep in mind,
% authors usually report the full width at half maximum (FWHM), which can
% be easily transformed to sigma.
% FWHM = sigma*sqrt(8*ln(2)) = sigma*2.354
target_duration = 314000; % 04:26 is the shortest scan duration; here in ms
scan_duration   = 400000;
TR              = 2.340;
% OJO: TR for KKI_1 EPI scans is 2340ms whereas for GU_1 it's 2530ms. TR
% should be in seconds, not miliseconds!

fil            = [0.01, 0.08]; % Filtering bandwidth
gsr            = 1;            % 1 = yes and 0 = no
wmcsfr         = 1;         % 1 = yes and 0 = no
WM_probseg     = 'anat/sub-01_space-MNI152NLin6Asym_res-2_label-WM_probseg.nii.gz';
CSF_probseg    = 'anat/sub-01_space-MNI152NLin6Asym_res-2_label-CSF_probseg.nii.gz';

%% START ----------------------------------------------------------------
% Note, in MATLAB, you can call a function without keyword assignments!

% Ensuring matching orientation
reorient_atlas(subjectdir, diratlas, atlas_filename, atlas);
atlas_filename = [diratlas 'reorient_fan246.nii.gz'];

%Trimming the Functional Scan to Duration 04:26 min
trimming(subjectdir, ...
         fsldir, ...
         filename, ...
         target_duration, ...
         scan_duration);

% Spatial Smoothing
spatial_smoothing(subjectdir, ...
                  fsldir, ...
                  smoothing_sigma);

filename = 'f0_smoothed';
% Temporal Filtering and Signal Regression
filtering_signalreg(subjectdir, ...
                    filename, ...
                    TR, ...
                    fil, ...
                    gsr, ...
                    wmcsfr, ...
                     WM_probseg, ...
                    CSF_probseg);

filename = 'f0_smoothed_fil_gsr_wmcsfr';
% Brain Parcellation and Z-scoring
functional_parcellation_zscore_AFNI(subjectdir, ...
                                    atlas_filename, ...
                                    filename, ...
                                    atlas);

%% ------------------------------------------------------------------------
% PARALLEL PROCESSING
% -------------------------------------------------------------------------
% %%%% Generates list of subjects to iterate over
% n_subj_st      = find(strcmp({dirdata.name}, 'sub-01')==1); % 1st subj
% n_subj_ed      = find(strcmp({dirdata.name}, 'sub-03')==1);
% % For n_subj_end: please put the largest subj number here.
% subj_ct=0;
% for s_row=n_subj_st:n_subj_ed
%     subj_ct=subj_ct+1;
%     subjs(subj_ct,:) = string([dirdata(s_row).folder, '/', char(dirdata(s_row).name)]);
% end
% Nsubjs=length(subjs);

% % Loop over participant directories
% for i = 1:length(subjs)
%     % The current subject's directory
%     current_subjectdir = subjs(i);
%     % Find the folder with functional scans in the participant's directory
%     func_dir   = fullfile(current_subjectdir, 'func');
%     % Get the files in the directory
%     func_files = dir(fullfile(func_dir, '*_task-rest_run-1_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz'));
%     % Apply the additional preprocessing functions in parallel
%     parfor x = 1:length(subjs)
%     % Spatial Smoothing
%     spatial_smoothing(subjectdir=current_subjectdir, ...
%                       filename=char(func_files.name), ...
%                       fsldir=fsldir, ...
%                       smoothing_sigma=smoothing_sigma);
%     end
% end

%% End --------------------------------------------------------------------
time = whatsthetime();
fprintf([time,' ', sub,'Finished the final preprocessing steps.\n'])
warning('on','all');