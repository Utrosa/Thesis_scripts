% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                        %
%                      Trimming the Functional Scans                     %
%                                                                        %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
function trimming(subjectdir, fsldir, filename, target_duration, ...
                  scan_duration)
% The durations should be in milliseconds.

input_file  = [subjectdir, 'func/', filename];
output_file = [subjectdir, 'f0_trimmed.nii.gz'];
k           = strfind(subjectdir,'sub');
sub         = [subjectdir(k:k+5),' - ']; clear k;
time        = whatsthetime();
fprintf([time, ' ', sub, 'Trimming the functional scan ... ']);

% Find the total number of timepoints
cmd = [fsldir, 'fslnvols ', input_file];

[~,timepoints] = system(cmd);
timepoints = str2double(timepoints);

% Calculate the number of the desired timepoints
new_duration = round((timepoints*target_duration)/scan_duration);

% Trim the scan
% Structure: fslroi <input> <output> <tmin> <tsize>
cmd = [fsldir, 'fslroi ', input_file, ' ', output_file, ...
       ' 0 ', num2str(new_duration)];
system(cmd);

fprintf('completed successfully.\n')
clear timepoints;