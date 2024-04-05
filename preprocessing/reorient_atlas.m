% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                        %
%                          Reorienting the Atlas                         %
%                                                                        %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
function reorient_atlas(subjectdir, diratlas, atlas_filename, atlas)

output_name = [diratlas 'reorient_' atlas '.nii.gz'];
k           = strfind(subjectdir,'sub');
sub         = [subjectdir(k:k+5),' - ']; clear k;
time        = whatsthetime();
fprintf([time, ' ', sub, 'Reorienting the atlas ... ']);

% Reorient along x dimension to standard: Left-to-Right
cmd = ['3dresample -orient LPI ' '-prefix ' output_name ...
       ' -input ' atlas_filename];
system(cmd);
fprintf('completed successfully.\n')