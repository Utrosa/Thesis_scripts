% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                        %
%                  Brain Parcellation: 3dROIstats in AFNI                %
%                                                                        %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
function functional_parcellation_zscore_AFNI(subjectdir, atlas_filename, ...
                                             filename, atlas)
% Extract the name of the subject that will be printed with every update.
k           = strfind(subjectdir,'sub');
sub         = [subjectdir(k:k+5),' - ']; 
clear k

% Specify the prerequisites 
cd(subjectdir);
time    = whatsthetime();

%% Start the Parcellation -------------------------------------------------
% Display current time, subject name, scan number, and preprocessing step 
% while the function is running
fprintf([time,' ',sub,': ','Functional Parcellation and Z-scoring ... '])
filename1 = [filename, '.nii.gz'];

% Compute the mean using only non_zero voxels -----------------------------
cmd_nzm = ['3dROIstats -mask ' atlas_filename ' -nomeanout -nzmean -quiet ' ...
            filename1 '  > ' filename '_seed_' atlas '.txt'];
system(cmd_nzm);

% Z-scoring ---------------------------------------------------------------
ts_parcel  = readmatrix([filename '_seed_' atlas '.txt'])';
ts_zscore  = zscore(ts_parcel,[],2);

% Compute the number of non_zero voxels: atlas ----------------------------
cmd_nzv_atlas = ['3dROIstats -mask ' atlas_filename ' -nomeanout -nzvoxels -quiet ' ...
                  atlas_filename '  > ' atlas '_voxelNum.txt'];
system(cmd_nzv_atlas);
voxel_per_roi1 = readmatrix([atlas '_voxelNum.txt']);

% Compute the number of non_zero voxels: filename -------------------------
% Capital V outputs voxel volume for non-zero voxels only.
cmd_nzv_data = ['3dROIstats -mask ' atlas_filename ' -nomeanout -nzvoxels -quiet ' ...
                 filename1 '[0]  > ' filename '_voxelNum_' atlas '.txt'];

system(cmd_nzv_data);
voxel_per_roi2    = readmatrix([filename '_voxelNum_' atlas '.txt']);
voxel_per_roi     = [voxel_per_roi1, voxel_per_roi2];

% Save the output ---------------------------------------------------------
save([filename '_on_parcel_' atlas '.mat'],'ts_parcel', ...
     'ts_zscore','voxel_per_roi');
fprintf('completed successfully.\n')