%% FINAL PREPROCESSING STEPS: Filtering & Signal Regression %%

% Set working directory
addpath("C:\Users\monik\Documents\fMRIprep_test\Flanker\derivatives\sub-08\func")

% Unzip the fMRI file
compressed_file = 'sub-08_task-flanker_run-2_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz';
unzipped_file = gunzip(compressed_file, '/Users/monik/Documents/fMRIprep_test/Flanker/derivatives/sub-08/func/');

% Load fMRI data
fmri_data = niftiread(unzipped_file{1});
%% 
% Define parameters
TR = 2; % Repetition time (in seconds)
nyquist_freq = 1 / (2 * TR); % Nyquist frequency
low_cutoff = 0.01; % Low cutoff frequency (in Hz)
high_cutoff = 0.08; % High cutoff frequency (in Hz)

% Perform FFT temporal filtering
n_volumes = size(fmri_volumes, 4);
n_voxels = numel(fmri_volumes(:, :, :, 1));

% Create frequency axis
freq_axis = (0:n_volumes-1) / (n_volumes * TR);

% Apply FFT to each voxel's time course
for voxel_idx = 1:n_voxels
    voxel_time_course = reshape(fmri_volumes(voxel_idx), [], n_volumes);
    voxel_fft = fft(voxel_time_course);
    
    % Apply bandpass filter
    voxel_fft(freq_axis < low_cutoff | freq_axis > high_cutoff, :) = 0;
    
    % Inverse FFT to get filtered time course
    filtered_time_course = real(ifft(voxel_fft));
    
    % Update voxel time course
    fmri_volumes(voxel_idx) = reshape(filtered_time_course, 1, 1, 1, n_volumes);
end

% Signal regression: global, white matter, and CSF signals regressed
% Here, you would implement the regression using the extracted signals
% For demonstration purposes, let's assume we have already extracted the
% global, white matter, and CSF signals and stored them in variables
% global_signal, white_matter_signal, and csf_signal respectively.

% Perform Z-score normalization on voxel time courses
mean_fmri_volumes = mean(fmri_volumes, 4);
std_fmri_volumes = std(fmri_volumes, [], 4);

zscored_fmri_volumes = (fmri_volumes - mean_fmri_volumes) ./ std_fmri_volumes;

% Now you have z-scored voxel time courses in zscored_fmri_volumes
% You can proceed with further analysis or save the data if needed

% For example, to save the z-scored fMRI data:
zscored_fmri_data = fmri_data;
zscored_fmri_data.img = zscored_fmri_volumes;
save_untouch_nii(zscored_fmri_data, 'zscored_fmri_data.nii');