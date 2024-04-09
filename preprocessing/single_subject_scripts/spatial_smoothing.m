% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                        %
%                       Spatial Smoothing: fslmaths                      %
%                                                                        %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
function spatial_smoothing(subjectdir, fsldir, smoothing_sigma)

%%% Define paths. Mark the subject and current time -----------------------
input_file  = [subjectdir, 'f0_trimmed.nii.gz'];
output_file = [subjectdir, 'f0_smoothed.nii.gz'];
k           = strfind(subjectdir,'sub');
sub         = [subjectdir(k:k+5),' - ']; clear k;
time        = whatsthetime();
fprintf([time, ' ', sub, 'Spatial smoothing ... ']);

%%% Construct the FSLmaths command ----------------------------------------
% Command structure: fslmaths input_image -s <smoothing sigma> output_image
fsl_cmd = sprintf('"%sfslmaths" "%s" -s %d "%s"', ...
                  fsldir, input_file, smoothing_sigma, output_file);

% Execute the FSL command
[status, cmdout] = system(fsl_cmd, '-echo');

% Check if command executed successfully
if status == 0
    disp('completed successfully.');
else
    error('Error:\n%s', cmdout);
end