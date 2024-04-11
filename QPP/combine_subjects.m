% ----------------------------------------------------------------------- %
%       COMIBINE PREPROCESSING OUTPUT: D0 MATRIX AS REQUIRED BY QPP       %
% ----------------------------------------------------------------------- %
clear, clc
% Get a list of subjects
subjects = {'sub-01', 'sub-02', 'sub-03'};

% Initialize a cell array
combined = cell(numel(subjects), 1);

% Loop through all participants
for i = 1:numel(subjects)

    % Define the directory for the current subject
    dirhead = pwd;
    subjectdir = [dirhead '\Input\' subjects{i}];

    % Find the file for the current subject
    filepath = fullfile(subjectdir, sprintf('%s_wmcsfr_gsr_fil_on_parcel_fan246.mat', subjects{i}));
    loadedfile = load(filepath); % ts_parcel or ts_zscore .. ?
    ts_parcel = loadedfile.ts_parcel;
    ts_parcel(1,:) = [];
    % Put the file in the matrix
    combined{i} = ts_parcel;
end

% Save
save('./Input/preprocessed_combined_test.mat', 'combined');