% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                        %
%             Brain Parcellation: FSLmeants for single ROI               %
%                                                                        %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
function functional_parcellation_FSL(subjectdir, atlas_filename, ...
                                     filename, atlas)
% Extract the name of the subject that will be printed with every update.
k           = strfind(subjectdir,'sub');
sub         = [subjectdir(k:k+5),' - ']; 
clear k

% Specify the prerequisites 
cd(subjectdir);
newStr  = split(filename,'_');
scan    = newStr{1};
time    = whatsthetime();

% Display current time, subject name, scan number, and preprocessing step 
% while the function is running
fprintf([time, ' ', sub, scan ': ', ...
         'Functional Parcellation ... '])

%% Start the Parcellation -------------------------------------------------
% fslmeants needs to know that the 4D image (input) is in the NifTi format
% and compressed. Without the appropriate extensions, files won't be read.
filename1 = [filename, '.nii.gz'];

% Resample mask matrix size to match input image size ---------------------
% cmd_resample = ['flirt -in ' filename1 ...
%                 ' -out resampled_f0.nii.gz ' ...
%                 ' -ref ' atlas_filename ' -applyxfm '];
% system(cmd_resample);

% Compute the mean using only non_zero voxels -----------------------------
cmd_nzmean = ['fslmeants' ' -i ' filename ...
              ' -o ' filename '_seed_' atlas '.txt' ...
              ' -m ' atlas_filename];
system(cmd_nzmean);

% Z-scoring
ts_parcel  = readmatrix([filename '_seed_' atlas '.txt'])';
ts_zscore  = zscore(ts_parcel,[],2);

% Compute the number of non_zero voxels: atlas ----------------------------
% Capital V outputs voxel volume for non-zero voxels only.
cmd_nzvoxel_atlas = ['fslstats ' atlas_filename ...
                     ' -k ' atlas_filename ' -V ' ...
                     '> ' atlas '_voxelNum.txt'];

system(cmd_nzvoxel_atlas);
voxel_per_roi1    = readmatrix([atlas '_voxelNum.txt']);

% Compute the number of non_zero voxels: filename -------------------------
% Capital V outputs voxel volume for non-zero voxels only.
cmd_nzvoxel_data  = ['fslstats ' filename1 ...
                     ' -k ' atlas_filename ' -V ' ... 
                     ' > ' filename '_voxelNum_' atlas '.txt'];

system(cmd_nzvoxel_data);
voxel_per_roi2    = readmatrix([filename '_voxelNum_' atlas '.txt']);
voxel_per_roi     = [voxel_per_roi1, voxel_per_roi2];

% Save the output ---------------------------------------------------------
save([filename '_on_parcel_' atlas '.mat'],'ts_parcel', ...
     'ts_zscore','voxel_per_roi');
fprintf('Completed functional parcellation for single ROI.\n')