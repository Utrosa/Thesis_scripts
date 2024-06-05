% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sig. results: which ROIs, confidence intervals, effect sizes, ...
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear; clc; close all;
%% Prepare the prerequisites
dirhead='...\Analysis\';
datadir=[dirhead, 'final_results\Static_FC\'];
% Q-values at FDR < 0.05
Qsig_matrix = load([datadir, 'Qvalues_sig_matrix.mat']);
Qsig_matrix = Qsig_matrix.q_values_thresholded;
% Correlation per subject
corr_matrix_ASD = load([datadir, 'subjects_correlation_matrix_ASD.mat']);
corr_matrix_TD  = load([datadir, 'subjects_correlation_matrix_TD.mat']);
corr_matrix_ASD = corr_matrix_ASD.correlationMatrices;
corr_matrix_TD  = corr_matrix_TD.correlationMatrices;
% Get a list of subjects: sublist_groupname.txt
group = {'ASD','TD'};
subjects_ASD = {};
subjects_TD = {};
for i = 1:length(group)
    subs_path = sprintf('%s\\data_static\\sublist_%s.txt',dirhead,group{i});
    subs_list = fopen(subs_path, 'r');
    subs = fscanf(subs_list, "%c");
    fclose(subs_list);
    subjects = regexp(subs, 'sub-\d{5}', 'match');
    if group{i} == "ASD"
        subjects_ASD = subjects;
    elseif group{i} == "TD"
        subjects_TD = subjects;
    end
end
clear subs; clear subs_list; clear subs_path;
% Get atlas parameters as in Xu et al. (2023)
AtlasTable = readtable([dirhead, 'data_static\resources\QPP_atlas.csv']);
Label='Label'; 
System='network7_shortname';
AtlasTable.Label=eval(['AtlasTable.' Label]); AtlasTable.System=eval(['AtlasTable.' System]);
AtlasTable = sortrows(AtlasTable,'Label','ascend');
[NetLB, net_index, ROI2Net]=unique(AtlasTable.System,'stable'); 
nnet=length(unique(ROI2Net)); iROI2NetC=cell(nnet,1); iNetL=zeros(nnet+1,1);
[a,iROI2Net]=sort(ROI2Net); 
for inet=1:nnet
    iROI2NetC{inet}=iROI2Net(a==inet);
    iNetL(inet+1)=iNetL(inet)+length(iROI2NetC{inet});
end
clear Label; clear a; clear System;
%% Find the atlas info about ROIs with significant static FC difference
% Find NaN & non-NaN values in Qsig_matrix
[nan_rows, nan_cols] = find(isnan(Qsig_matrix));
[non_nan_rows, non_nan_cols] = find(~isnan(Qsig_matrix));
num_non_nan = length(non_nan_rows);
Qvalue_at_non_nan_loc = zeros(num_non_nan, size(Qsig_matrix,3));

% These you found from the Qsig_Matrix were reorder based on their location
% in Yeo's 7 network (indices correspond to iROI2Net). It is therefore very
% important to also reorder the atlas accordingly.
AtlasTable = AtlasTable(iROI2Net, :);
%% Loop through significant connections to obtain Q-values
for k = 1:num_non_nan
    i = non_nan_rows(k);
    j = non_nan_cols(k);
    Qvalue_at_non_nan_loc(k, :) = Qsig_matrix(i, j, :);
end

% Loop through all sig. differences to find ROI indices
pairs_idx = cell(length(non_nan_rows), 1);
for k = 1:length(non_nan_rows)
    pairs_idx{k} = {non_nan_rows(k), non_nan_cols(k)};
end
% Find the brain area names for the ROI pairs
pairs_ROIs = cell(length(non_nan_rows), 1); % Short anatomical label
pairs_desc = cell(length(non_nan_rows), 1); % Long anatomical label
pairs_net7 = cell(length(non_nan_rows), 1); % Yeo's 7 Network name
pairs_net17 = cell(length(non_nan_rows), 1); % Yeo's 17 Network name
for k = 1:length(non_nan_rows)
    row_id = pairs_idx{k}{1};
    col_id = pairs_idx{k}{2};
    % Access the relevant variables in the Atlas spreadsheet
    pairs_ROIs{k} = [AtlasTable.region{row_id}, '; ', AtlasTable.region{col_id}];
    pairs_desc{k} = [AtlasTable.area_description{row_id}, '; ', AtlasTable.area_description{col_id}];
    pairs_net7{k} = [AtlasTable.network7_shortname{row_id}, '; ', AtlasTable.network7_shortname{col_id}];
    pairs_net17{k} = [AtlasTable.network17_shortname{row_id}, '; ', AtlasTable.network17_shortname{col_id}];
end
clear col_id; clear row_id;
%% Find the significant static FC values in subjects' matrices
num_non_nan = length(non_nan_rows);
ASD_corr_at_non_nan_loc = zeros(num_non_nan, size(corr_matrix_ASD,3));
TD_corr_at_non_nan_loc  = zeros(num_non_nan, size(corr_matrix_TD, 3));
% Loop through significant connections
for k = 1:num_non_nan
    i = non_nan_rows(k);
    j = non_nan_cols(k);
    ASD_corr_at_non_nan_loc(k, :) = corr_matrix_ASD(i, j, :);
    TD_corr_at_non_nan_loc(k, :) = corr_matrix_TD(i, j, :);
end
%% Calculate descriptives for significant results
% Mean static FC (correlation coefficient) for each pair of regions
ASD_mean = mean(ASD_corr_at_non_nan_loc,2);
TD_mean  = mean(TD_corr_at_non_nan_loc,2);

% Standard deviation for each pair of regions
ASD_std = std(ASD_corr_at_non_nan_loc,0,2);
TD_std  = std(TD_corr_at_non_nan_loc,0,2);

% Effect size: Cohen's d
effect_size = [];
confidence_interval1 = [];
confidence_interval2 = [];
for i = 1:length(ASD_mean)
    es = meanEffectSize(ASD_corr_at_non_nan_loc(i, :), ...
                        TD_corr_at_non_nan_loc(i, :),Effect="cohen");
    effect_size = [effect_size, es.Effect];
    confidence_interval1 = [confidence_interval1, es.ConfidenceIntervals(1)];
    confidence_interval2 = [confidence_interval2, es.ConfidenceIntervals(2)];
end
%% Structure the output in a table
ASD_pairs = array2table(ASD_corr_at_non_nan_loc,'VariableNames',subjects_ASD(:));
TD_pairs  = array2table(TD_corr_at_non_nan_loc, 'VariableNames',subjects_TD);

ASD_mean = array2table(ASD_mean,'VariableNames',{'mean_staticFC_ASD'});
TD_mean  = array2table(TD_mean, 'VariableNames',{'mean_staticFC_TD'});

ASD_std = array2table(ASD_std, 'VariableNames',{'std_staticFC_ASD'});
TD_std  = array2table(TD_std, 'VariableNames',{'std_staticFC_TD'});

Q_values  = array2table(Qvalue_at_non_nan_loc, 'VariableNames',{'q_value'});

effect_size = array2table(effect_size', 'VariableNames',{'cohen_d'});
confidence_interval1 = array2table(confidence_interval1', 'VariableNames',{'conf_int1'});
confidence_interval2 = array2table(confidence_interval2', 'VariableNames',{'conf_int2'});

% Create tables which you can save as excel spreadsheets
ROI_info = table(pairs_idx, pairs_ROIs, pairs_desc, pairs_net7, pairs_net17);
result_subjects_ASD = [ROI_info, ASD_pairs];
result_subjects_TD  = [ROI_info, TD_pairs];

result_group = [ROI_info, ASD_mean, TD_mean, ASD_std, TD_std, Q_values, ...
                effect_size, confidence_interval1, confidence_interval2];
%% SAVE THE RESULTS
writetable(result_subjects_ASD,[datadir,'staticFC_subASD.xlsx'])
writetable(result_subjects_TD,[datadir,'staticFC_subTD.xlsx'])
writetable(result_group,[datadir,'staticFC_group.xlsx'])