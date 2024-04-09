% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                        %
%                          Reorienting the Atlas                         %
%                                                                        %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
function reorient_atlas(diratlas, atlas_filename, atlas)

% Define the name of the reoriented atlas (output).
output_name = [diratlas 'reorient_' atlas '.nii.gz'];
time        = whatsthetime();
fprintf([time, ' Reorienting the atlas ... ']);

% Reorient along x dimension to standard: Left-to-Right
cmd = ['3dresample -orient LPI ' '-prefix ' output_name ...
       ' -input ' atlas_filename];
system(cmd);
fprintf('completed successfully.\n')