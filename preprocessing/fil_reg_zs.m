% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
%      Filtering, Signal Regression, Brain Parcellation, Z-scoring        %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
clear; close all; clc

% Setup -------------------------------------------------------------------
addpath('C:\Program Files\MATLAB\R2022b');
addpath('C:\Users\monik\NIfTI_20140122');

dirhead     = pwd;
dirdata     = dir([dirhead '/test-data']);
%n_subj_st   = find(strcmp({dirdata.name}, 'subject008')==1);

subjectdir     = [dirhead '/test-data/subject008/'];
diratlas       = [dirhead '/preprocessing_Xu_Smith_2023/resources'];
atlas          = 'fan246'; %Brainnetome atlas 2016 with 246 parcels
atlas_filename = [diratlas '/Fan2016parcel_Yeo/BN_Atlas_246_2mm.nii.gz']; 

% Input the smoothed fMRIPrep preprocessed data
filename = 'f01_smoothed';

% Set the parameters: TR (depends on ABIDE collection), fill (filtering 
% bandwidth), gsr (global signal regression; 1 = yes and 0 = no), wmcsfr 
% (white matter and cerebrospinal fluid regression; 1 = yes and 0 = no).
TR = 2;
fil = [0.01, 0.08];
gsr = 1;
wmcsfr = 1;

% Which subject?
k = strfind(subjectdir,'subject');
sub = [subjectdir(k:end),' - ']; 
clear k

warning('off','all');

% Temporal Filtering of Functional Data -----------------------------------
time = whatsthetime();
fprintf([time,' ',sub,'Temporal filtering ... '])

% Load data
cd (subjectdir)
nii = load_untouch_nii([filename,'.nii.gz']);
f = nii.img;
[X,Y,Z,T] = size(f);

% Getting the dimensions of the functional scan
f = reshape(f,[X*Y*Z,T]); 
f = f';

% Reshaping the scan into a 2D matrix of time over voxels
tc = f;

% Renaming the functional timeseries as tc
SF = 1/TR';

% Calculating the sampling frequency using the TR
x = size(tc,1);

% x is the number of timepoints
tc_mirror = tc(end:-1:1,:);

% This the mirror of the timeseries
tc_extended = [tc_mirror; tc; tc_mirror];

% Extend the timecourse in a periodical way
tc_fft = fft(tc_extended);

% Getting the fast fourier transform
if fil(1) > 0
    high_cutoff = round(size(tc_extended,1)*fil(1)/SF);
    tc_fft(1:high_cutoff,:) = 0;
end

% Conducting the high pass filtering
if fil(2) > 0
    low_cutoff = round(size(tc_extended,1)*fil(2)/SF);
    tc_fft(low_cutoff:end,:) = 0;
end

% Conducting the low pass filtering
tc_filtered = ifft(tc_fft);

% Inverse fast fourier transform
temp = 2*real(tc_filtered);
tc_filtered = temp(x+1:2*x,:);

% Getting the final filtered functional scan
f = tc_filtered';
f = reshape(f,[X,Y,Z,T]);
nii.img = f;
save_untouch_nii(nii,[filename,'_fil.nii.gz']);

% Just to cut down on file size, I'm multiplying the functional
% scan by the brain mask again
filename = [filename,'_fil'];
fprintf('Done\n')
clear f1 tc_filtered temp tc_fft low_cutoff high_cutoff
clear tc_extended tc_mirror x

% Global signal regression ------------------------------------------------
%
time = whatsthetime();
fprintf([time,' ', sub, 'Global signal regression ... '])
if gsr == 1
    nii = load_untouch_nii([filename,'.nii.gz']);
    f = nii.img;
    % Loading the functional scan into Matlab
    f_gsr = f;
    % Predefining the global signal regressed functional scan
    f_mean = zeros(1,size(f,4));
    % Predefining what will me the mean of the functional scan
    for t = 1:size(f,4)
        temp = f(:,:,:,t);
        temp_sum = sum(temp(:));
        temp_length = length(find(temp(:)));
        temp_mean = temp_sum/temp_length;
        f_mean(t) = temp_mean;
    end
    % Taking the mean of the functional timeseries. This is the
    % functional timecourse.
    f_mean_zsc = zscore(f_mean);
    % Z-scoring the mean
    for x = 1:size(f,1)
        for y = 1:size(f,2)
            for z = 1:size(f,3)
                v = double(squeeze(f(x,y,z,:)));
                beta = (f_mean_zsc * f_mean_zsc') \ ...
                    (f_mean_zsc * v);
                f_gsr(x,y,z,:) = v' - f_mean_zsc * beta;
            end
        end
    end
    % Regressing the global signal from the functional
    % timeseries
    f = f_gsr; 
    filename = [filename,'_gsr'];
    nii.img = f;
    save_untouch_nii(nii,[filename,'.nii.gz']);
    % Saving the global signal regressed functional scan
    fprintf('Done\n')
else
    fprintf('Skipped\n')
end

% White Matter and CSF signal regression ----------------------------------
time = whatsthetime();
fprintf([time, ' ', sub, 'WM/CSF signal regression ... '])
if wmcsfr == 1
    if gsr == 0
        nii = load_untouch_nii([filename,'.nii.gz']);
        f = nii.img;
        % Loading the functional scan
    end
    wm_mask = load_nii('sub-08_space-MNI152NLin2009cAsym_res-2_label-WM_probseg.nii.gz');
    wm_mask = wm_mask.img;
    % Loading white matter mask
    csf_mask = ...
        load_nii('sub-08_space-MNI152NLin2009cAsym_res-2_label-CSF_probseg.nii.gz');
    csf_mask = csf_mask.img;
    % Loading the CSF mask
    wmcsf_mask = double(logical(wm_mask + csf_mask));
    % Creating a mask of the white matter and CSF together and
    % making sure that it is a binary matrix
    f_wmcsf = double(f);
    % Predefining the matrix that will hold just the timeseries
    % from the white matter and CSF 
    for t = 1:size(f,4)
        f_wmcsf(:,:,:,t) = f_wmcsf(:,:,:,t) .* wmcsf_mask;
    end
    % Removing all the gray matter voxels from f_wmcsf
    f_wmcsf_r = f_wmcsf;
    % Predefining the matrix that will be the regressed signal
    % from the white matter and CSF
    f_wmcsf_mean = zeros(1,size(f,4));
    % Predefining what will me the mean of the functional scan
    for t = 1:size(f,4)
        temp = f_wmcsf(:,:,:,t);
        temp_sum = sum(temp(:));
        temp_length = length(find(temp(:)));
        temp_mean = temp_sum/temp_length;
        f_wmcsf_mean(t) = temp_mean;
    end
    % Taking the mean of the wmcsf signal. This is the
    % functional timecourse of the wmcsf.
    f_wmcsf_mean_zsc = zscore(f_wmcsf_mean);
    % Z-scoring the mean
    for x = 1:size(f,1)
        for y = 1:size(f,2)
            for z = 1:size(f,3)
                v = double(squeeze(f(x,y,z,:)));
                beta = (f_wmcsf_mean_zsc*f_wmcsf_mean_zsc') ...
                    \ (f_wmcsf_mean_zsc * v);
                f_wmcsf_r(x,y,z,:) = ...
                    v' - f_wmcsf_mean_zsc * beta;
            end
        end
    end
    % Regressing the white matter and CSF signal from the
    % functional timeseries 
    f = f_wmcsf_r;
    nii.img = f;
    save_untouch_nii(nii,[filename,'_wmcsfr.nii.gz']);
    % Saving the white matter/CSF regressef functiona scan
    filename = [filename,'_wmcsfr'];
    fprintf('Done\n')
else
    fprintf('Skipped\n')
end

% Brain Parcellation ------------------------------------------------------
if gsr==1
    ext='gsr';
    filename1=[filename '_' ext];
    functional_parcellation_zscore(subjectdir, atlas_filename, filename1, atlas);
end

if wmcsfr==1
    ext='wmcsfr';
    filename1=[filename '_' ext];
    functional_parcellation_zscore(subjectdir, atlas_filename, filename1, atlas)
end
% Z-scoring all voxels in functional scan ______________________________ %
%
time = whatsthetime();
fprintf([time,' ', sub, 'Z-scoring ... '])
if gsr ~= 1
    nii = load_untouch_nii([filename,'.nii.gz']);
    f = nii.img;
    % Loading the functional scan in case it wasn't already
    % loaded in the previous step
end
mask = load_nii([atlasdir,'/atlas_t1_brain_mask.nii.gz']);
mask = mask.img;
% Loading the T1 mask
[x,y,z,~] = size(f);
% Assigning variables to the dimensions of the scan
for j = 1:x
    for k = 1:y
        for l = 1:z
            if mask(j,k,l) == 1
                voxel = f(j,k,l,:);
                % Working on one brain voxel at a time
                voxel_zsc = zscore(voxel);
                % Z-scoring the voxel signal
                f(j,k,l,:) = voxel_zsc;
                % Adding the zscored value to the zscored
                % functional scan matrix
            end
        end
    end
end
nii.img = f;
save_untouch_nii(nii,[filename,'_zsc.nii.gz']);
filename = [filename,'_zsc'];
clear nii x y z j k l
fprintf('Done\n')
% Saving the newly z-scored functional scan
% ___________________________________________________________ %

time = whatsthetime();
fprintf([time,' ', sub,'Finished the final preprocessing steps\n'])
warning('on','all');