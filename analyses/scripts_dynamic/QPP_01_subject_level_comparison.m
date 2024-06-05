% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This script compares the similarity between subjects' QPPs
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear; clc; close all;
%% Prepare the prerequisites
% Set the data directories
datadir = '.../Analysis/';
GU_folder_path   = 'data_dynamic/QPP_GU/Output/SbjQPP/';
KKI_folder_path  = 'data_dynamic/QPP_KKI/Output/SbjQPP/';
OHSU_folder_path = 'data_dynamic/QPP_OHSU/Output/SbjQPP/';
% Get a list of subjects: sublist_groupname.txt
group = {'ASD','TD'};
subjects_ASD = {};
subjects_TD = {};
for i = 1:length(group)
    subs_path = sprintf('%sdata_static/sublist_%s.txt',datadir,group{i});
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
% Get file names for each collection
GU_files = dir(fullfile([datadir, GU_folder_path], '*.mat'));
KKI_files = dir(fullfile([datadir, KKI_folder_path], '*.mat'));
OHSU_files = dir(fullfile([datadir, OHSU_folder_path], '*.mat'));
clear OHSU_folder_path; clear KKI_folder_path; clear GU_folder_path;
% Join all file names into one struct
all_files = vertcat(GU_files, KKI_files, OHSU_files);
%% Load all QPP output files per subject
all_data = struct();
for i = 1:length(all_files)
    file_path = fullfile(all_files(i).folder, all_files(i).name);
    data = load(file_path);
    all_data(i).name = all_files(i).name;
    all_data(i).data = data;
end
clear all_files; clear data; clear i; clear file_path;
%% Prepare the Sliding Correlation Vectors (Cs & Crs)
% Remove trailing zeros from correlation vectors
for i = 1:length(all_data)
    % Select the sliding correlation vectors for each participants
    cs_vector = all_data(i).data.Cs;
    cs_vector_reg = all_data(i).data.Crs{1,1};
    % Find the last non-zero element
    last_nonzero_idx1 = find(cs_vector, 1, 'last');
    last_nonzero_idx2 = find(cs_vector_reg, 1, 'last');
    % Trim the vector to remove trailing zeros
    trimmed_cs_vector = cs_vector(1:last_nonzero_idx1);
    trimmed_cs_vector_reg = cs_vector_reg(1:last_nonzero_idx2);
    % Update the original structure with the trimmed vector
    all_data(i).data.Cs = trimmed_cs_vector;
    all_data(i).data.Crs{1,1} = trimmed_cs_vector_reg;
end
clear cs_vector; clear last_nonzero_index1; clear trimmed_cs_vector;
clear cs_vector_reg; clear last_nonzero_index2; clear trimmed_cs_vector_reg;
%% Calculate Basic QPP Metrics: QPP strength and frequency per subject
% QPP strength is defined as the mean peak height in the sliding
% correlation vectors (Cs), wheras QPP frequency is defined as the rate of 
% peak occurrence.

%%% Function definition: estimate strength and frequency for subjects' QPPs
%%%% Function input should be an array: data.Cs or data.Crs{1,1}
function [mean_peak_height, peaks, locations, rate] = QPP_metrics(data)
    cs_vector = abs(data);
    [peaks, locations] = findpeaks(cs_vector,"MinPeakHeight",0.2);
    mean_peak_height = mean(peaks);
    rate = length(peaks) / length(cs_vector);
end

%%% Compute basic QPP metrics for all subjects' correlation vectors
num_files = length(all_data);
%%%% Initialize empty arrays
QPP_metrics_ASD = struct();
QPP_metrics_TD  = struct();
ASD_Cs = struct(); ASD_Crs = struct();
TD_Cs = struct(); TD_Crs = struct();
%%%% Counter
asd_count = 0;
td_count = 0;
%%%% Loop through all participants
for i = 1:num_files
    data = all_data(i).data;
    cs_vector = all_data(i).data.Cs;
    cs_vector_reg = all_data(i).data.Crs{1,1};
    [peaks_cs, ~,~, rate_cs] = QPP_metrics(data.Cs);
    [peaks_reg,~,~,rate_reg] = QPP_metrics(data.Crs{1,1});
    % Compute QPP strength and frequency for ASD
    if contains(all_data(i).name, 'ASD')
        asd_count = asd_count + 1;
        QPP_metrics_ASD(asd_count).strengthCs  = peaks_cs;
        QPP_metrics_ASD(asd_count).strengthCrs = peaks_reg;
        QPP_metrics_ASD(asd_count).frequencyCs = rate_cs;
        QPP_metrics_ASD(asd_count).frequencyCrs = rate_reg;
        ASD_Cs(asd_count).Cs   = cs_vector;
        ASD_Crs(asd_count).Crs = cs_vector_reg;
    % Compute QPP strength and frequency for TD
    elseif contains(all_data(i).name, 'TD')
        td_count = td_count + 1;
        QPP_metrics_TD(td_count).strengthCs = peaks_cs;
        QPP_metrics_TD(td_count).frequencyCs = rate_cs;
        QPP_metrics_TD(td_count).strengthCrs = peaks_reg;
        QPP_metrics_TD(td_count).frequencyCrs = rate_reg;
        TD_Cs(td_count).Cs   = cs_vector;
        TD_Crs(td_count).Crs = cs_vector_reg;
    end
end
clear asd_count; clear td_count; clear data; clear num_files; clear i;
clear cs_vector; clear cs_vector_reg; clear rate_cs; clear rate_reg;
clear peaks_cs; clear peaks_reg;

%%% Add subject IDs to the QPP Metrics struct
%%%% ASD Participants
for i = 1:numel(QPP_metrics_ASD)
    QPP_metrics_ASD(i).subjectID = subjects_ASD{i};
end
%%%% Control Participants
for i = 1:numel(QPP_metrics_TD)
    QPP_metrics_TD(i).subjectID = subjects_TD{i};
end
%%% Save vectors of the basic QPP metrics for each subject
% save([datadir, 'final_results/Dynamic_FC/cs_vector_ASD.mat'],'ASD_Cs');
% save([datadir, 'final_results/Dynamic_FC/cs_vector_TD.mat'],'TD_Cs');
% save([datadir, 'final_results/Dynamic_FC/crs_vector_ASD.mat'],'ASD_Crs');
% save([datadir, 'final_results/Dynamic_FC/crs_vector_TD.mat'],'TD_Crs');
% save([datadir, 'final_results/Dynamic_FC/QPP_metrics_ASD.mat'],'QPP_metrics_ASD');
% save([datadir, 'final_results/Dynamic_FC/QPP_metrics_TD.mat'],'QPP_metrics_TD');

%%% Save also an Excel Spreadsheet of the Subjects' Metrics
QPP_metrics_ASD_table = struct2table(QPP_metrics_ASD);
writetable(QPP_metrics_ASD_table,[datadir,'final_results/Dynamic_FC/QPP_metrics_ASD.xlsx']);

QPP_metrics_TD_table = struct2table(QPP_metrics_TD);
writetable(QPP_metrics_TD_table,[datadir, 'final_results/Dynamic_FC/QPP_metrics_TD.xlsx']);
%% Plotting 
QPP_metrics_ASD = struct2table(QPP_metrics_ASD);
QPP_metrics_TD  = struct2table(QPP_metrics_TD);
%%% Replace NaNs with 0 (otherwise no plot possible & Nan == 0 here)
idx1 = ismissing(QPP_metrics_ASD(:,{'strengthCrs','frequencyCrs'}));
QPP_metrics_ASD{:,{'strengthCrs','frequencyCrs'}}(idx1) = 0;
idx2 = ismissing(QPP_metrics_TD(:,{'strengthCrs','frequencyCrs'}));
QPP_metrics_TD{:,{'strengthCrs','frequencyCrs'}}(idx2) = 0;

%%% QPP Strength histograms per subject group
figure;
histogram(QPP_metrics_ASD.strengthCs,'FaceColor','blue','FaceAlpha',0.4,'EdgeColor','black');hold on;
histogram(QPP_metrics_TD.strengthCs,'FaceColor','red','FaceAlpha',0.4,'EdgeColor','black');
xlabel('Peak Height','FontName','Calibri','FontSize',12);
ylabel('Frequency','FontName','Calibri','FontSize',12);
title('QPP Strength','FontName','Calibri','FontSize',14);
legend({'ASD', 'TD'},'FontName','Calibri','FontSize',10);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w'); box off; grid off; hold off;
f = gcf;
exportgraphics(f,[datadir,'\final_results\Dynamic_FC\figs\QPP_strength_group.png'],'Resolution',300)
savefig(f,[datadir,'\final_results\Dynamic_FC\figs\QPP_strength_group.fig']);

%%% QPP Frequency histograms per subject group
figure;
histogram(QPP_metrics_ASD.frequencyCs,'FaceColor','blue','FaceAlpha',0.4,'EdgeColor','black');hold on;
histogram(QPP_metrics_TD.frequencyCs,'FaceColor','red','FaceAlpha',0.4,'EdgeColor','black');
xlabel('Peak Height','FontName','Calibri','FontSize',12);
ylabel('Frequency','FontName','Calibri','FontSize',12);
title('QPP Frequency','FontName','Calibri','FontSize',14);
legend({'ASD', 'TD'},'FontName','Calibri','FontSize',10);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w'); box off; grid off; hold off;
f = gcf;
exportgraphics(f,[datadir,'\final_results\Dynamic_FC\figs\QPP_frequency_group.png'],'Resolution',300)
savefig(f,[datadir,'\final_results\Dynamic_FC\figs\QPP_frequency_group.fig']);

%%% QPP strength vectors: before and after regression
figure;
%%%% ASD participants
subplot(1,2,1);
histogram(QPP_metrics_ASD.strengthCs,'FaceColor','blue','EdgeColor','black');hold on;
histogram(QPP_metrics_ASD.strengthCrs,'FaceColor','red','EdgeColor','black');
xlabel('Peak Height','FontName','Calibri','FontSize',12);
ylabel('Frequency','FontName','Calibri','FontSize',12);
title('ASD Participants','FontName','Calibri','FontSize',14);
legend({'Before regression', 'After regression'},'FontName','Calibri','FontSize',10);
set(gca,'FontName','Calibri','FontSize',10); hold off;
set(gcf,'Color','w'); box off; grid off;

%%%% Control participants
subplot(1,2,2);
histogram(QPP_metrics_TD.strengthCs,'FaceColor','blue','EdgeColor','black');hold on;
histogram(QPP_metrics_TD.strengthCrs,'FaceColor','red','EdgeColor','black');
xlabel('Peak Height','FontName','Calibri','FontSize',12);
ylabel('Frequency','FontName','Calibri','FontSize',12);
title('Control Participants','FontName','Calibri','FontSize',14);
legend({'Before regression', 'After regression'},'FontName','Calibri','FontSize',10);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w'); box off; grid off; hold off;

f = gcf;
exportgraphics(f,[datadir,'\final_results\Dynamic_FC\figs\QPP_strength_before_after_reg.png'],'Resolution',300)
savefig(f,[datadir,'\final_results\Dynamic_FC\figs\QPP_strength_before_after_reg.fig']);

%%% QPP frequency vectors per subjects
figure;
%%%% ASD participants
subplot(1,2,1);
histogram(QPP_metrics_ASD.frequencyCs,'FaceColor','blue','EdgeColor','black');hold on;
histogram(QPP_metrics_ASD.frequencyCrs,'FaceColor','red','EdgeColor','black');
xlabel('Peak Rate','FontName','Calibri','FontSize',12);
ylabel('Frequency','FontName','Calibri','FontSize',12);
title('ASD Participants','FontName','Calibri','FontSize',14);
legend({'Before regression', 'After regression'},'FontName','Calibri','FontSize',10);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w'); box off; grid off; hold off;

%%%% Control participants
subplot(1,2,2);
histogram(QPP_metrics_TD.frequencyCs,'FaceColor','blue','EdgeColor', 'black');hold on;
histogram(QPP_metrics_TD.frequencyCrs,'FaceColor','red','EdgeColor', 'black');
xlabel('Peak Rate','FontName','Calibri','FontSize',12);
ylabel('Frequency','FontName','Calibri','FontSize',12);
title('Control Participants','FontName','Calibri','FontSize',14);
legend({'Before regression', 'After regression'},'FontName','Calibri','FontSize',10);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w'); box off; grid off; hold off;

f = gcf;
exportgraphics(f,[datadir,'\final_results\Dynamic_FC\figs\QPP_frequency_before_after_reg.png'],'Resolution',300)
savefig(f,[datadir,'\final_results\Dynamic_FC\figs\QPP_frequency_before_after_reg.fig']);

%% QPP frequency & mean for ASD
figure;
%%% ASD participants
subplot(1,2,1);
means1=[mean(QPP_metrics_ASD.strengthCs),mean(QPP_metrics_ASD.strengthCrs);...
        mean(QPP_metrics_ASD.frequencyCs),mean(QPP_metrics_ASD.frequencyCrs)];
labels = ["Strength" "Frequency"];
bar(labels,means1);colororder("glow");
legend({'Before regression', 'After regression'},'FontName','Calibri','FontSize',10);
ylabel('Mean correlation/frequency of peaks','FontName','Calibri','FontSize',12);
title('ASD Participants','FontName','Calibri','FontSize',14);
set(gca,'FontName','Calibri','FontSize',10);set(gcf,'Color','w'); box off; grid off;

%%% Control participants
subplot(1,2,2);
means2=[mean(QPP_metrics_TD.strengthCs),mean(QPP_metrics_TD.strengthCrs);...
        mean(QPP_metrics_TD.frequencyCs),mean(QPP_metrics_TD.frequencyCrs)];
bar(labels,means2);colororder("glow");
ylabel('Mean correlation/frequency of peaks','FontName','Calibri','FontSize',12);
title('Control Participants','FontName','Calibri','FontSize',14);
legend({'Before regression', 'After regression'},'FontName','Calibri','FontSize',10);
set(gca,'FontName','Calibri','FontSize',10);set(gcf,'Color','w'); box off; grid off;

f = gcf;
exportgraphics(f,[datadir,'\final_results\Dynamic_FC\figs\QPP_group_before_after_reg.png'],'Resolution',300)
savefig(f,[datadir,'\final_results\Dynamic_FC\figs\QPP_group_before_after_reg.fig']);
%% Perform two-sample t-tests 
%%%% Group difference in QPP metrics (ASD vs TD)
[~, p_stcs, ci_stcs, stats_stcs] = ttest2(QPP_metrics_ASD.strengthCs,...
                                          QPP_metrics_TD.strengthCs);
[~, p_frcs, ci_frcs, stats_frcs] = ttest2(QPP_metrics_ASD.frequencyCs,...
                                          QPP_metrics_TD.frequencyCs);
%%%% Differece in QPP metrics before and after regression
[~, p_stcrsA, ci_stcrsA, stats_stcrsA] = ttest2(QPP_metrics_ASD.strengthCs,...
                                      QPP_metrics_ASD.strengthCrs);
[~, p_stcrsT, ci_stcrsT, stats_stcrsT] = ttest2(QPP_metrics_TD.strengthCs,...
                                      QPP_metrics_TD.strengthCrs);
% Display the QPP stats 
fprintf('QPP STRENGTH\n');
fprintf('Mean strength ASD: %f\n', mean(QPP_metrics_ASD.strengthCs));
fprintf('Mean strength TD: %f\n', mean(QPP_metrics_TD.strengthCs));
fprintf('Standard deviation ASD: %f\n', std(QPP_metrics_ASD.strengthCs));
fprintf('Standard deviation TD: %f\n', std(QPP_metrics_TD.strengthCs));

fprintf('Two-sample t-test results:\n');
fprintf('p-value: %.4f\n', p_stcs);
fprintf('Confidence Interval: [%.4f, %.4f]\n', ci_stcs(1), ci_stcs(2));
fprintf('Test Statistic: %.4f\n', stats_stcs.tstat);
clear p_stcs; clear ci_stcs; clear stats_stcs;

% Display the QPP frequency results
fprintf('\nQPP FREQUENCY\n');
fprintf('Mean frequency ASD: %f\n', mean(QPP_metrics_ASD.frequencyCs));
fprintf('Mean frequency TD: %f\n', mean(QPP_metrics_TD.frequencyCs));
fprintf('Standard deviation ASD: %f\n', std(QPP_metrics_ASD.frequencyCs));
fprintf('Standard deviation TD: %f\n', std(QPP_metrics_TD.frequencyCs));

fprintf('Two-sample t-test results:\n');
fprintf('p-value: %.4f\n', p_frcs);
fprintf('Confidence Interval: [%.4f, %.4f]\n', ci_frcs(1), ci_frcs(2));
fprintf('Test Statistic: %.4f\n', stats_frcs.tstat);
clear p_frcs; clear ci_frcs; clear stats_frcs;

% Display differences in QPP basic metrics before and after regression
fprintf('\nQPP REGRESSION\n');

fprintf('STRENGTH ASD\n');
fprintf('Two-sample t-test results:\n');
fprintf('p-value: %.4f\n', p_stcrsA);
fprintf('Confidence Interval: [%.4f, %.4f]\n', ci_stcrsA(1), ci_stcrsA(2));
fprintf('Test Statistic: %.4f\n', stats_stcrsA.tstat);
clear p_stcrs; clear ci_stcrs; clear stats_stcrs;

fprintf('STRENGTH TD\n');
fprintf('Two-sample t-test results:\n');
fprintf('p-value: %.4f\n', p_stcrsT);
fprintf('Confidence Interval: [%.4f, %.4f]\n', ci_stcrsT(1), ci_stcrsT(2));
fprintf('Test Statistic: %.4f\n', stats_stcrsT.tstat);
clear p_stcrsT; clear ci_stcrsT; clear stats_stcrsT; clear i;