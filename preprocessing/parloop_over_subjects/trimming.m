% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                        %
%                      Trimming the Functional Scans                     %
%                                                                        %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
function trimming(subjectdir, fsldir, filename, target_duration, ...
                  scan_duration)

%%% Prepare the input file: fMRIPrep-preprocessed functional scans.
input_file  = [subjectdir, '/func/', filename];
output_file = [subjectdir, '_trimmed.nii.gz']; 
% The output will be saved in the derivatives folder.

%%% Get information about subject and current time
k           = strfind(subjectdir,'sub');
sub         = [subjectdir(k:k+5),' - ']; clear k;
time        = whatsthetime();
fprintf([time, ' ', sub, 'Trimming the functional scan ... ']);

%%% Find the total number of timepoints
cmd = [fsldir, 'fslnvols ', input_file];
[~,timepoints] = system(cmd);
timepoints = str2double(timepoints);

%%% Calculate the number of the desired timepoints
new_duration = round((timepoints*target_duration)/scan_duration);

%%% Trim the scan
% FSL command structure: fslroi <input> <output> <tmin> <tsize>
cmd = [fsldir, 'fslroi ', input_file, ' ', output_file, ...
       ' 0 ', num2str(new_duration)];
system(cmd);

fprintf('completed successfully.\n')
clear timepoints;