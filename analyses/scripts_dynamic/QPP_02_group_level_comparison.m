% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This script compares the similarity between group correlation vectors
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear; clc; close all;

%%% SET THE PARAMETERS
% This script requires consistency in Scan vs Subject QPP output files, as
% well as the type of QPP we're looking at: QPP, phase-adjusted, or
% phase-revered. The combination: subject & QPPas is not optimal.
cs_type = 'scan'; % Less interpolation (concatenation per site)
% cs_type = 'subject'; % More interpolation (for every subject)
% qpp = 'QPPs'; % QPPs(1) is the QPP, QPPs(2) is the phase-reversed QPP
qpp = 'QPPas'; % phase-adjusted QPP

% Set common time bins. The sites collected functional data with 
% different TR (repetition time), so we now need to get a common time 
% axis for the QPP outputs. time = TR * timepoints (bins)
common_QPP_timepoints = 80;
common_QPP_time = linspace(0, 40, common_QPP_timepoints);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Setup
%%% Set the paths
datadir = '.../Analysis/';
if cs_type=="scan"
    GU_folder_path = 'data_dynamic/QPP_GU/Output/ScanQPP/';  % Concatenated
    KKI_folder_path = 'data_dynamic/QPP_KKI/Output/ScanQPP/';
    OHSU_folder_path = 'data_dynamic/QPP_OHSU/Output/ScanQPP/';
elseif cs_type=="subject"
    GU_folder_path = 'data_dynamic/QPP_GU/Output/SbjQPP/';
    KKI_folder_path = 'data_dynamic/QPP_KKI/Output/SbjQPP/';
    OHSU_folder_path = 'data_dynamic/QPP_OHSU/Output/SbjQPP/';
end

%%% Get file names
GU_files = dir(fullfile([datadir, GU_folder_path], '*.mat'));
KKI_files = dir(fullfile([datadir, KKI_folder_path], '*.mat'));
OHSU_files = dir(fullfile([datadir, OHSU_folder_path], '*.mat'));
all_files = vertcat(GU_files, KKI_files, OHSU_files);
clear GU_folder_path; clear KKI_folder_path, clear OHSU_folder_path;
clear KKI_files; clear GU_files; clear OHSU_files;

%%% Join files per diagnostic group
ASD_files = struct();
TD_files = struct();
asd_count=0;
td_count=0;
for i = 1:length(all_files)
    file_path = fullfile(all_files(i).folder, all_files(i).name);
    data = load(file_path);

    if contains(all_files(i).name, 'ASD')
        asd_count = asd_count + 1;
        ASD_files(asd_count).data = data;
        ASD_files(asd_count).file_name = all_files(i).name;
    elseif contains(all_files(i).name, 'TD')
        td_count = td_count + 1;
        TD_files(td_count).data = data;
        TD_files(td_count).file_name = all_files(i).name;
    end
end
clear asd_count; clear td_count; clear data; clear i; clear file_path;
clear all_files;

%%% Save
save([datadir, 'final_results/Dynamic_FC/ASD_files.mat'], 'ASD_files');
save([datadir, 'final_results/Dynamic_FC/TD_files.mat'], 'TD_files');
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Function Definition
%%% Define functions to interpolate QPPs/as matrices to a common X axis
function interpolated_matrices = interpolate_QPPs(...
    input_files,common_QPP_timepoints,common_QPP_time)

    % Initialize empty cell array for the interpolated matrices
    interpolated_matrices = cell(length(input_files), 1);
    
    % Loop through the original QPP matrices
    for i = 1:length(input_files)
        [num_regions, num_timepoints] = size(input_files(i).data.QPPs{1});
        original_time = linspace(0, 40, num_timepoints);
        
        % Interpolate each row (region) to the common time vector
        interpolated_matrix = zeros(num_regions, common_QPP_timepoints);
        for region = 1:num_regions
            interpolated_matrix(region, :) = interp1( ...
            original_time, input_files(i).data.QPPs{1}(region, :), ...
            common_QPP_time, 'makima');
        end
        interpolated_matrices{i} = interpolated_matrix;
    end
end

function interpolated_matrices = interpolate_QPPas(...
    input_files,common_QPP_timepoints,common_QPP_time)
    % Initialize empty cell array for the interpolated matrices
    interpolated_matrices = cell(length(input_files), 1);
    % Loop through the original QPP matrices
    for i = 1:length(input_files)
        [num_regions, num_timepoints] = size(input_files(i).data.QPPas{1}{1});
        original_time = linspace(0, 40, num_timepoints);
        
        % Interpolate each row (region) to the common time vector
        interpolated_matrix = zeros(num_regions, common_QPP_timepoints);
        for region = 1:num_regions
            interpolated_matrix(region, :) = interp1( ...
            original_time, input_files(i).data.QPPas{1}{1}(region, :), ...
            common_QPP_time, 'makima');
        end
        interpolated_matrices{i} = interpolated_matrix;
    end
end
%% Create Group QPP
%%% Interpolate QPPs/as matrices
if qpp=="QPPs"
    interpolated_ASD = interpolate_QPPs(ASD_files,common_QPP_timepoints,common_QPP_time);
    interpolated_TD  = interpolate_QPPs(TD_files,common_QPP_timepoints,common_QPP_time);
elseif qpp=="QPPas"
    interpolated_ASD = interpolate_QPPas(ASD_files,common_QPP_timepoints,common_QPP_time);
    interpolated_TD  = interpolate_QPPas(TD_files,common_QPP_timepoints,common_QPP_time);
end

%%% Filter out empty arrays (this occurs for cs_type=="subject")
non_empty_ASD = interpolated_ASD(~cellfun('isempty', interpolated_ASD));
non_empty_TD  = interpolated_TD(~cellfun('isempty',  interpolated_TD));
clear interpolated_ASD; clear interpolated_TD;

%%% Concatenate non-empty arrays along the third dimension
%%%% ASD
if ~isempty(non_empty_ASD)
    all_ASD_matrices = cat(3, non_empty_ASD{:});
end
%%%% Control
if ~isempty(non_empty_TD)
    all_TD_matrices = cat(3, non_empty_TD{:});
end
clear non_empty_TD;clear non_empty_ASD;clear common_QPP_timepoints;

%%% Obtain the group QPPs/as matrix by averaging
average_ASD_matrix = mean(all_ASD_matrices, 3);
average_TD_matrix  = mean(all_TD_matrices, 3);

%%% Get 2D QPP (whole-brain)
groupASD_QPP = mean(average_ASD_matrix(ASD_files(1).data.iROI2Net,:)', 2);
groupTD_QPP  = mean(average_TD_matrix(TD_files(1).data.iROI2Net,:)', 2);

%%% Plot the group QPPs/as timecourses (over time in seconds)
figure;
plot(common_QPP_time, groupASD_QPP,'b','LineWidth', 1.5); hold on; 
plot(common_QPP_time, groupTD_QPP,'r','LineWidth', 1.5);box off; grid off;
title('Group QPP','FontName','Calibri','FontSize',14);
xlabel('Time (s)','FontName','Calibri','FontSize', 12);
ylabel('Mean Signal Amplitude','FontName','Calibri','FontSize', 12);
legend({'ASD', 'TD'}, 'FontName', 'Calibri', 'FontSize', 10);
set(gca, 'FontName','Calibri','FontSize', 10);set(gcf,'Color','w');
f = gcf;
exportgraphics(f,[datadir,'/final_results/Dynamic_FC/figs/group_',qpp,'_',cs_type,'_averaged.png'],'Resolution',300)
savefig(f,[datadir,'/final_results/Dynamic_FC/figs/group_',qpp,'_',cs_type,'_averaged.fig']);
clear f;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Group Sliding Correlation Vector Before (Cs) and Asfter (Crs) Regression
%%% Functions to interpolate subject and scan vectors
function interpolated_subjects = interpolate_Sbj(input_files,common_timepoints,common_time)
    % Initialize empty cell array for the interpolated vectors
    interpolated_subjects = cell(length(input_files), 1);
    % Loop through the original vectors
    for i = 1:length(input_files)
        [num_subjects, num_timepoints] = size(input_files(i).Cs);
        original_time = linspace(0, 280, num_timepoints);
        % Interpolate each row (subject) to the common time vector
        interpolated_subject = zeros(num_subjects, common_timepoints);
        for sub = 1:num_subjects
            interpolated_subject(sub, :) = interp1( ...
            original_time,input_files(i).Cs(sub, :),common_time,'makima');
        end
        interpolated_subjects{i} = interpolated_subject;
    end
end
function interpolated_subjects = interpolate_SbjReg(input_files,common_timepoints,common_time)
    % Initialize empty cell array for the interpolated vectors
    interpolated_subjects = cell(length(input_files), 1);
    % Loop through the original vectors
    for i = 1:length(input_files)
        [num_subjects, num_timepoints] = size(input_files(i).Crs);
        original_time = linspace(0, 280, num_timepoints);
        % Interpolate each row (subject) to the common time vector
        interpolated_subject = zeros(num_subjects, common_timepoints);
        for sub = 1:num_subjects
            interpolated_subject(sub, :) = interp1( ...
            original_time,input_files(i).Crs(sub, :),common_time,'makima');
        end
        interpolated_subjects{i} = interpolated_subject;
    end
end
function interpolated_scans = interpolate_scan(input_files,common_timepoints,common_time)
    interpolated_scans = cell(length(input_files), 1);
    % Loop through the original scan vectors
    for i = 1:length(input_files)
        % Remove trailing zeros & update the file
        cs_vector = input_files(i).data.Cs;
        last_nonzero_index = find(cs_vector, 1, 'last');
        trimmed_cs_vector  = cs_vector(1:last_nonzero_index);
        input_files(i).data.Cs = trimmed_cs_vector;
        % Interpolate
        [~, num_timepoints] = size(input_files(i).data.Cs);
        original_time = linspace(0, 6240, num_timepoints);
        interpolated_scans{i} = interp1( ...
        original_time,input_files(i).data.Cs(1, :),common_time,'makima');
    end
end
function interpolated_scans = interpolate_scanReg(input_files,common_timepoints,common_time)
    interpolated_scans = cell(length(input_files), 1);
    % Loop through the original scan vectors
    for i = 1:length(input_files)
        % Remove trailing zeros & update the file
        crs_vector = input_files(i).data.Crs{1};
        last_nonzero_index = find(crs_vector, 1, 'last');
        trimmed_crs_vector  = crs_vector(1:last_nonzero_index);
        input_files(i).data.Crs{1} = trimmed_crs_vector;
        % Interpolate
        [~, num_timepoints] = size(input_files(i).data.Crs{1});
        original_time = linspace(0, 6240, num_timepoints);
        interpolated_scans{i} = interp1( ...
        original_time,input_files(i).data.Crs{1}(1, :),common_time,'makima');
    end
end

%%% Load or compute the group Cs & Crs (before and after regression)
if cs_type == "subject"
    % Load the averaged cs_vectors which we created using the script 
    % QPP_01_subject_level_comparison.m (Cs & Crs)
    data_ASD = load(fullfile(datadir,'final_results/Dynamic_FC/cs_vector_ASD.mat'));
    data_ASD = data_ASD.ASD_Cs;
    data_ASD_reg = load(fullfile(datadir,'final_results/Dynamic_FC/crs_vector_ASD.mat'));
    data_ASD_reg = data_ASD_reg.ASD_Crs;
    data_TD  = load(fullfile(datadir,'final_results/Dynamic_FC/cs_vector_TD.mat'));
    data_TD  = data_TD.TD_Cs;
    data_TD_reg = load(fullfile(datadir,'final_results/Dynamic_FC/crs_vector_TD.mat'));
    data_TD_reg = data_TD_reg.TD_Crs;
    % Subject Cs vectors are shorter (around 5 min per subject)
    common_timepoints = 113; % The min no. of time bins across collections
    common_time = linspace(0, 280, common_timepoints);
    % Interpolate files
    interpolated_vec_ASD = interpolate_Sbj(data_ASD,common_timepoints,common_time);
    interpolated_vec_TD  = interpolate_Sbj(data_TD,common_timepoints,common_time);
    interpolated_vec_ASD_reg = interpolate_SbjReg(data_ASD_reg,common_timepoints,common_time);
    interpolated_vec_TD_reg  = interpolate_SbjReg(data_TD_reg,common_timepoints,common_time);

elseif cs_type == "scan"
    data_ASD = ASD_files;
    data_TD  = TD_files;
    % Cs vectors from Scan output are longer than vectors from Sbj output!
    common_timepoints = 3120; % The min no. of time bins across collections
    common_time = linspace(0, 6240, common_timepoints);
    % Interpolate files
    interpolated_vec_ASD = interpolate_scan(data_ASD,common_timepoints,common_time);
    interpolated_vec_TD  = interpolate_scan(data_TD,common_timepoints,common_time);
    interpolated_vec_ASD_reg = interpolate_scanReg(data_ASD,common_timepoints,common_time);
    interpolated_vec_TD_reg  = interpolate_scanReg(data_TD,common_timepoints,common_time);
end

%%% Concatenate the interpolated Cs & Crs vectors into a 3D array
all_ASD_vectors = cat(3, interpolated_vec_ASD{:});
all_TD_vectors  = cat(3, interpolated_vec_TD{:});
all_ASD_vectors_reg = cat(3, interpolated_vec_ASD_reg{:});
all_TD_vectors_reg  = cat(3, interpolated_vec_TD_reg{:});

%%% Compute the average Cs & Crs
average_ASD_vector = mean(all_ASD_vectors, 3);
average_TD_vector = mean(all_TD_vectors, 3);
average_ASD_vector_reg = mean(all_ASD_vectors_reg, 3);
average_TD_vector_reg = mean(all_TD_vectors_reg, 3);
clear common_timepoints;

%%% Plot Group Sliding Correlation Vector
figure;
plot(common_time,average_ASD_vector,'b','LineWidth',1);
hold on; box off; grid off;
plot(common_time,average_TD_vector,'r','LineWidth',1);
xlabel('Time (s)', 'FontName', 'Calibri','FontSize',12);%xlim([5900 6200])
ylabel('Mean Correlation', 'FontName', 'Calibri','FontSize',12);
title('Group Sliding Correlation','FontName', 'Calibri','FontSize',14);
legend({'ASD', 'TD'}, 'FontName', 'Calibri', 'FontSize', 10);
set(gca,'FontName','Calibri','FontSize',10); set(gcf,'Color','w');hold off;
f = gcf;
exportgraphics(f,[datadir,'/final_results/Dynamic_FC/figs/',qpp,'_',cs_type,'Cs_average.png'],'Resolution',300)
savefig(f,[datadir,'/final_results/Dynamic_FC/figs/',qpp,'_',cs_type,'Cs_average.fig']);
clear f;

%% Descriptive Statistics of the Group Correlation Vectors (Cs & Crs)
%%% Compute strength of QPP over the GROUP correlation vector (average)
[peaks_ASD, locations_ASD] = findpeaks(abs(average_ASD_vector),"MinPeakHeight",0.2);
[peaks_TD, locations_TD] = findpeaks(abs(average_TD_vector),"MinPeakHeight",0.2);
[peaks_ASD_reg, locations_ASD_reg] = findpeaks(abs(average_ASD_vector_reg),"MinPeakHeight",0.2);
[peaks_TD_reg, locations_TD_reg] = findpeaks(abs(average_TD_vector_reg),"MinPeakHeight",0.2);

mean_peak_height_ASD = mean(peaks_ASD);
mean_peak_height_TD  = mean(peaks_TD);
mean_peak_height_ASD_reg = mean(peaks_ASD_reg);
mean_peak_height_TD_reg  = mean(peaks_TD_reg);

no_peaks_ASD = length(peaks_ASD);
no_peaks_TD  = length(peaks_TD);
no_peaks_ASD_reg = length(peaks_ASD_reg);
no_peaks_TD_reg  = length(peaks_TD_reg);

fprintf('Group QPP strength in ASD: %.4f\n',mean_peak_height_ASD);
fprintf('Group QPP strength in TD:  %.4f\n',mean_peak_height_TD);

fprintf('Standard deviation ASD: %f\n', std(peaks_ASD));
fprintf('Standard deviation TD: %f\n', std(peaks_TD))

fprintf('\nNumber of peaks in group QPP in ASD: %.4f\n',no_peaks_ASD);
fprintf('Number of peaks in group QPP in TD:  %.4f\n',no_peaks_TD);

%%% Compute frequency of QPP over the GROUP correlation vector (average)
rate_ASD = no_peaks_ASD / length(average_ASD_vector);
rate_TD  = no_peaks_TD / length(average_TD_vector);
rate_ASD_reg = no_peaks_ASD_reg / length(average_ASD_vector_reg);
rate_TD_reg  = no_peaks_TD_reg / length(average_TD_vector_reg);

fprintf('\nFrequency of group QPP in ASD: %.4f\n',rate_ASD);
fprintf('Frequency of group QPP in TD:  %.4f\n',rate_TD);

%%% Compute the mean square difference in group QPP strength
squaredDifferences = (groupASD_QPP - groupTD_QPP).^2;
MSD = mean(squaredDifferences);
max_diff = max(squaredDifferences);
fprintf('\nGROUP QPP');
fprintf('\nMean Square Difference: %.4f\n',MSD);
fprintf('Max Square Difference: %.4f\n',max_diff);
fprintf('Square root of Mean Square Difference: %.4f\n',sqrt(MSD));

% Compute the correlation between the signals
correlationCoefficient = corr(groupASD_QPP, groupTD_QPP);
fprintf('Correlation Coefficient: %.4f\n', correlationCoefficient);

%%% Compute the mean square difference in group sliding correlation vectors
squaredDifferences = (average_ASD_vector - average_TD_vector).^2;
MSD = mean(squaredDifferences);
max_diff = max(squaredDifferences);
fprintf('\nGROUP Cs');
fprintf('\nMean Square Difference: %.4f\n',MSD);
fprintf('Max Square Difference: %.4f\n',max_diff);
fprintf('Square root of Mean Square Difference: %.4f\n',sqrt(MSD));

%% Inferential Statistics on Group Outputs
%%% Two-sample t test on the group QPPas (phase-adjusted > QPPs)
[~, p_qpp, ci_qpp, stats_qpp] = ttest2(groupASD_QPP,groupTD_QPP);
fprintf('\nGROUP QPP - Two-sample t-test results:\n');
fprintf('p-value: %.4f\n', p_qpp);
fprintf('Confidence Interval: [%.4f, %.4f]\n', ci_qpp(1), ci_qpp(2));
fprintf('Test Statistic: %.4f\n', stats_qpp.tstat);
clear p_qpp; clear ci_qpp; clear stats_qpp;

%%% Two-sample t test on the group sliding correlation vectors (Cs)
[~, p_cs, ci_cs, stats_cs] = ttest2(average_ASD_vector,average_TD_vector);
fprintf('\nGROUP SLIDING CORRELATION - Two-sample t-test results:\n');
fprintf('p-value: %.4f\n', p_cs);
fprintf('Confidence Interval: [%.4f, %.4f]\n', ci_cs(1), ci_cs(2));
fprintf('Test Statistic: %.4f\n', stats_cs.tstat);
clear p_cs; clear ci_cs; clear stats_cs;

%%% Two-sample t test on the basic metrics of Cs
[~, p_stcs, ci_stcs, stats_stcs] = ttest2(peaks_ASD,peaks_TD);
fprintf('\nBasic Metrics - GROUP Cs STRENGTH\n');
fprintf('Two-sample t-test results:\n');
fprintf('p-value: %.4f\n', p_stcs);
fprintf('Confidence Interval: [%.4f, %.4f]\n', ci_stcs(1), ci_stcs(2));
fprintf('Test Statistic: %.4f\n', stats_stcs.tstat);
%% Plotting the difference in peak magnitude (QPP strength)
%%% Sort the peaks in descending order
pks_sorted_ASD = findpeaks(abs(average_ASD_vector),"MinPeakHeight",0.2,'SortStr','descend');
pks_sorted_TD  = findpeaks(abs(average_TD_vector),"MinPeakHeight",0.2,'SortStr','descend');
pks_sorted_ASD = pks_sorted_ASD(:);
pks_sorted_TD  = pks_sorted_TD(:);

%%% Create a common x-axis index
numPeaks = max(no_peaks_ASD, no_peaks_TD);
x = 1:numPeaks;

%%% Prepare data for plotting, padding with NaNs where necessary
pks_ASD_padded = [pks_sorted_ASD; NaN(numPeaks - no_peaks_ASD, 1)];
pks_TD_padded  = [pks_sorted_TD; NaN(numPeaks - no_peaks_TD, 1)];
figure;
plot(x, pks_ASD_padded,'b','LineWidth', 1.5);hold on;box off;grid off;
plot(x, pks_TD_padded,'r', 'LineWidth', 1.5);
xlabel('Peak Index','FontSize',12);
ylabel('Peak Height (correlation)','FontSize',12);
title('Sorted Peaks','FontSize',14);
legend({'ASD', 'TD'}, 'FontName', 'Calibri', 'FontSize', 10);
set(gca, 'FontName', 'Calibri', 'FontSize', 10);
set(gcf, 'Color', 'w'); hold off;
f = gcf;
exportgraphics(f,[datadir,'/final_results/Dynamic_FC/figs/',qpp,'_',cs_type,'Cs_peaks_sorted.png'],'Resolution',300)
savefig(f,[datadir,'/final_results/Dynamic_FC/figs/',qpp,'_',cs_type,'Cs_peaks_sorted.fig']);

%%% Sliding correlation vectors: before and after regression
figure;
%%%% ASD participants
subplot(1,2,1);
histogram(average_ASD_vector,'FaceColor','blue','EdgeColor','black',FaceAlpha=0.4);hold on;
histogram(average_ASD_vector_reg,'FaceColor','red','EdgeColor','black',FaceAlpha=0.4);
xlabel('Sliding Correlation','FontName','Calibri','FontSize',12);
ylabel('Frequency','FontName','Calibri','FontSize',12);
title('ASD Participants','FontName','Calibri','FontSize',14);
legend({'Before QPP Removal', 'After QPP Removal'},'FontName','Calibri','FontSize',10);
set(gca,'FontName','Calibri','FontSize',10); hold off;
set(gcf,'Color','w'); box off; grid off;

%%%% Control participants
subplot(1,2,2);
histogram(average_TD_vector,'FaceColor','blue','EdgeColor','black',FaceAlpha=0.4);hold on;
histogram(average_TD_vector_reg,'FaceColor','red','EdgeColor','black',FaceAlpha=0.4);
xlabel('Sliding Correlation','FontName','Calibri','FontSize',12);
ylabel('Frequency','FontName','Calibri','FontSize',12);
title('Control Participants','FontName','Calibri','FontSize',14);
legend({'Before QPP removal', 'After QPP removal'},'FontName','Calibri','FontSize',10);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w'); box off; grid off; hold off;

f = gcf;
exportgraphics(f,[datadir,'\final_results\Dynamic_FC\figs\sliding_corr_before_after_reg.png'],'Resolution',300)
savefig(f,[datadir,'\final_results\Dynamic_FC\figs\sliding_corr_strength_before_after_reg.fig']);

%% Compute descriptive and inferential statistics of QPP metrics per scan
if cs_type=="scan"
    
    % Initialize arrays to store all peaks and locations for each group
    peaks_ASD_all = [];
    peaks_TD_all = [];
    rate_ASD_all = [];
    rate_TD_all = [];
    
    % Loop through each group
    groups = {'ASD', 'TD'};
    
    for g = 1:length(groups)
        group = groups{g};
        
        % Extract interpolated vector based on group
        interpolated_vec = eval(['interpolated_vec_' group]);
        
        % Loop through each scan
        peaks_group = [];
        rate_group = zeros(1,3);
        for i = 1:3
            % Find peaks and concatenate
            [peaks_scan, ~] = findpeaks(abs(interpolated_vec{i}), 'MinPeakHeight', 0.2);
            peaks_group =  [peaks_group; peaks_scan(:)];

            % Calculate frequency for each scan (don't concatenate)
            rate_scan = length(peaks_scan) / numel(interpolated_vec{i});
            rate_group(i) = rate_scan;
        end
        
        % Store all peaks and frequencies for the ASD & TD groups
        eval(['peaks_' group '_all = peaks_group;']);
        eval(['rate_' group '_all  = rate_group;']);
    end
    
    % Perform t-tests between ASD and TD for strength and frequency
    [~, p_strength]  = ttest2(peaks_ASD_all, peaks_TD_all);
    [~, p_frequency] = ttest2(rate_ASD_all,  rate_TD_all);
    
    % Print results for differences across dataset collections (n=3)
    fprintf('\nBasic Metrics - SCAN (SITE) Cs STRENGTH\n');
    fprintf('T-test for strength between ASD and TD: %.4f\n', p_strength);
    fprintf('T-test for frequency between ASD and TD: %.4f\n', p_frequency);
end