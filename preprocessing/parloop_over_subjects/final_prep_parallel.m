% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
%                 Smoothing, Filtering, Signal Regression,                %
%                    Brain Parcellation, and Z-scoring                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Note, this script will only work on a Linux system (or WSL for Windows).
clear; close all; clc
fprintf('\nBeginning the additional preprocessing pipeline\n\n');

%% SETUP ENVIRONMENT AND MAIN DIRECTORY -----------------------------------
% Add the path to Nifti toolbox, which is necessary to read nii.gz files.
% Add the path to the main directory, where the Atlas and ABIDE data 
% folders are located.
addpath("/USER/MATLAB/NIfTI_20140122/", ...
        "/USER/preprocessing/")

%%% Set the FSL environment and output filetype to '.nii.gz'
fsldir = '/usr/local/fsl/';
setenv('FSLDIR','/usr/local/fsl');
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

%%% Set the main directory
dirhead = pwd;
dirdata = dir([dirhead '/ABIDE_test/derivatives']); 
% ABIDE_test is the folder with test data. The folder has to be in BIDS.

%%% Find the atlas
diratlas = [dirhead '/Atlas/'];
atlas    = 'fan246'; %Brainnetome atlas 2016 with 246 parcels
% atlas_filename = [diratlas 'BN_Atlas_246_2mm.nii.gz'];
% reorient_atlas(diratlas, atlas_filename, atlas);
atlas_filename = [diratlas 'reorient_fan246.nii.gz'];
% To ensure a matching orientation between atlas and functional images run
% the reorient_atlas() function. If this has already been done, you can 
% skip this step and set the value of the atlas_filename variable directly 
% to the name of the reoriented atlas.

warning('off','all');
% There are a few warnings that MATLAB outputs depending on the type of the
% computer this function is being run on. As far as my understanding goes,
% they are irrelevant, hence the function is temporarily turning off the
% warnings, so that the updates being printed as the code runs look clean
% (see the Anzar Abbas' preprocessing pipeline from 2016).

%% SET THE PARAMETERS -----------------------------------------------------
%%% Trimming parameters
target_duration = 314000; % 04:26 is the shortest scan duration; here in ms
scan_duration   = 400000; % As reported by authors of KKI_1 (ABIDE-II).
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

%% SETUP PARTICIPANT DIRECTORIES AND FILENAMES ----------------------------
%%% Generates a subject list to iterate over (see Nan Xu's preprocessing
%%% pipeline from 2023).
n_subj_st      = find(strcmp({dirdata.name}, 'sub-01')==1); % 1st subj
n_subj_ed      = find(strcmp({dirdata.name}, 'sub-03')==1);
% For n_subj_end: please put the largest subj number here.
subjs = strings(n_subj_ed - n_subj_st + 1, 1);
subj_ct = 0;
for s_row = n_subj_st:n_subj_ed
    subj_ct = subj_ct+1;
    subjs(subj_ct,:) = string([dirdata(s_row).folder, '/', ...
                       char(dirdata(s_row).name)]);
end
Nsubjs = length(subjs);

%%% Loop over participant directories
filenames    = cell(length(subjs), 1);
WM_probsegs  = cell(length(subjs), 1);
CSF_probsegs = cell(length(subjs), 1);

for i = 1:length(subjs)
    % The current subject's directory
    current_subjectdir = subjs{i};
    % Find the folders with functional and anatomical scans
    func_dir = fullfile(current_subjectdir, 'func');
    anat_dir = fullfile(current_subjectdir, 'anat');
    % Get the preprocessed BOLD filenames in the directory
    preproc_bold = dir(fullfile(func_dir, '*_task-rest_run-1_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz'));
    % Get the white matter and cerebrospinal fluid masks
    WM_probseg = dir(fullfile(anat_dir, '*_space-MNI152NLin6Asym_res-2_label-WM_probseg.nii.gz'));
    CSF_probseg = dir(fullfile(anat_dir, '*_space-MNI152NLin6Asym_res-2_label-CSF_probseg.nii.gz'));
    % Store filenames in the preallocated cell array
    filenames{i} = preproc_bold.name;
    WM_probsegs{i} = WM_probseg.name;
    CSF_probsegs{i} = CSF_probseg.name;
end

%%% Convert to string array
filenames = string(filenames);
WM_probsegs = string(WM_probsegs);
CSF_probsegs = string(CSF_probsegs);

%% START ------------------------------------------------------------------
%%% Call the additional preprocessing functions in parallel
parfor subj_ct = 1:length(subjs)
    subjectdir = char(subjs(subj_ct));
    filename_subj = char(filenames(subj_ct,:));
    WM_probseg_subj = char(WM_probsegs(subj_ct,:));
    CSF_probseg_subj = char(CSF_probsegs(subj_ct,:));

    k = strfind(subjectdir,'sub');
    sub = [subjectdir(k:k+5)];
    % This is simply helping with the text that will be displayed as the
    % function is running. 'sub' is a string with the name of the subject 
    % that will be printed with every update (see the Anzar Abbas'
    % preprocessing pipeline from 2016).
    
    %Trimming the functional scan to specified duration
    trimming(subjectdir, ...
             fsldir, ...
             filename_subj, ...
             target_duration, ...
             scan_duration);

    % Spatial smoothing
    spatial_smoothing(subjectdir, ...
                      fsldir, ...
                      smoothing_sigma);

    % Temporal filtering and signal regression
    filtering_signalreg(subjectdir, ...
                        TR, ...
                        fil, ...
                        gsr, ...
                        wmcsfr, ...
                        WM_probseg_subj, ...
                        CSF_probseg_subj);

    % Brain parcellation with specified atlas and z-scoring
    functional_parcellation_zscore_AFNI(subjectdir, ...
                                        atlas_filename, ...
                                        atlas);
end

%% End --------------------------------------------------------------------
time = whatsthetime();
fprintf([time,' Finished the final preprocessing steps.\n'])
warning('on','all');