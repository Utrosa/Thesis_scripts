% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                        %
%     MATLAB script to perform smoothing of fMRI data using FSLmaths     %
%                                                                        %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Define paths ------------------------------------------------------------
dirhead = pwd;
input_file = [dirhead, '/subject008/f01.nii.gz'];
output_file = [dirhead, '/subject008/f01_smoothed_WSL.nii.gz'];
subject_dir = 'subject008/f01.nii.gz';

% Set the environment variables
fprintf('\nBeginning the spatial smoothing\n\n');
fsldir = '/usr/local/fsl/share/fsl/bin/';
setenv('FSLDIR','/usr/local/fsl');
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

% Add path to NifTi toolbox
addpath("/mnt/c/Users/monik/Documents/MATLAB/NIfTI_20140122/")

% Set the smoothing parameter: sigma --------------------------------------
% CAREFUL, authors usually report the full width at half maximum (FWHM).
% FWHM = sigma*sqrt(8*ln(2)) = sigma*2.354
smoothing_sigma = 6/sqrt(8*log(2));

% Who is the current participant? -----------------------------------------
k = strfind(subject_dir,'subject');  
sub = [subject_dir(k:k+9),' - ']; 
clear k

% What's the time now? ----------------------------------------------------
time = whatsthetime();

% Tell me who's data was smoothed when ------------------------------------
fprintf([time, ' ', sub, ' Spatial smoothing ... ']);

% Construct FSLmaths command ----------------------------------------------
% fslmaths input_image -s <smoothing sigma> output_image
fsl_cmd = sprintf('"%sfslmaths" "%s" -s %d "%s"', ...
                  fsldir, input_file, smoothing_sigma, output_file);

% Execute the command in WSL ----------------------------------------------
[status, cmdout] = system(fsl_cmd, '-echo');

% Check if command executed successfully
if status == 0
    disp('DONE with smoothing.');
else
    error('Error:\n%s', cmdout);
end
% ______________________________________________________________________ %