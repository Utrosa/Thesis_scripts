% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SAVE SPACE: COMBINE MATRICES
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear; clc;
dirhead='...\Analysis\';
datadir=[dirhead, '\final_results\Static_FC\'];
%% Load Matrices and Atlas Parameters
% fc_matrix_1 = load([datadir, 'averaged_correlation_matrix_ASD.mat']);
% fc_matrix_2 = load([datadir, 'averaged_correlation_matrix_TD.mat']);
% fc_matrix_1 = fc_matrix_1.averageCorrelationMatrix;
% fc_matrix_2 = fc_matrix_2.averageCorrelationMatrix;

% fc_matrix_1 = load([datadir, 'averaged_correlation_matrix_reg_ASD.mat']);
% fc_matrix_2 = load([datadir, 'averaged_correlation_matrix_reg_TD.mat']);
% fc_matrix_1 = fc_matrix_1.FC_REG_ASD;
% fc_matrix_2 = fc_matrix_2.FC_REG_TD;

% fc_matrix_1 = load([datadir, 'averaged_correlation_matrix_ASD.mat']);
% fc_matrix_2 = load([datadir, 'averaged_correlation_matrix_reg_ASD.mat']);
% fc_matrix_1 = fc_matrix_1.averageCorrelationMatrix;
% fc_matrix_2 = fc_matrix_2.FC_REG_ASD;

% fc_matrix_1 = load([datadir, 'Pvalues_matrix.mat']);
% fc_matrix_2 = load([datadir, 'Pvalues_reg_matrix.mat']);
% fc_matrix_1 = fc_matrix_1.p_values;
% fc_matrix_2 = fc_matrix_2.p_values;

% fc_matrix_1 = load([datadir, 'Qvalues_matrix.mat']);
% fc_matrix_2 = load([datadir, 'Qvalues_reg_matrix.mat']);
% fc_matrix_1 = fc_matrix_1.q_values;
% fc_matrix_2 = fc_matrix_2.q_values;

fc_matrix_1 = load([datadir, 'Qvalues_sig_matrix.mat']);
fc_matrix_2 = load([datadir, 'Qvalues_reg_sig_matrix.mat']);
fc_matrix_1 = fc_matrix_1.q_values_thresholded;
fc_matrix_2 = fc_matrix_2.q_values_thresholded;

% Get atlas paramaters as in Xu et al. (2023)
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
clear Label; clear System; clear a;
%% Combine Matrices
num_rois = size(fc_matrix_1,1);
combined_matrix = zeros(num_rois,num_rois);
for i = 1:num_rois
    combined_matrix(i+1:end,i) = fc_matrix_1(i+1:end,i);
    combined_matrix(i,i+1:end) = fc_matrix_2(i,i+1:end);
end
%% Plot the Combined Matrices
%%% Adjust the scale: for correlation [-1 1] and for p-values [0 1].

%%% Figure setting for matrices
% imagesc(combined_matrix, [-1 1]);
% colorbar;colormap(jet);plotNets(ROI2Net, NetLB, 10, 1);
% hold on;
% diag_x = 1:size(combined_matrix, 1);
% diag_y = 1:size(combined_matrix, 2);
% plot(diag_x, diag_y, 'white', 'LineWidth', 1); % 'k' for black color
% set(gca, 'FontName', 'Calibri', 'FontSize', 14);
% hold off;

%%% Figure settings for Q-values matrices with thresholded values (NaNs)
figure;
significance_level = 0.05; % FDR < 0.05
imagesc(combined_matrix,[0 significance_level]);
custom_colormap = [0.3 0.3 0.4; jet(256)]; % Add gray to the colormap
colormap(custom_colormap); plotNets(ROI2Net, NetLB, 10, 1);colorbar;
% title('Q-values (FDR < 0.05) After Regression','FontName','Calibri','FontSize',14);
hold on;
diag_x = 1:size(combined_matrix, 1);
diag_y = 1:size(combined_matrix, 2);
plot(diag_x, diag_y, 'white', 'LineWidth', 1); % 'k' for black color
set(gca, 'FontName', 'Calibri', 'FontSize', 14);
hold off;

% figure;
% significance_level = 0.2; % r > 0.1
% imagesc(combined_matrix,[0 significance_level]);
% custom_colormap = [0.5 0.5 0.5; jet(256)]; % Add gray to the colormap
% colormap(custom_colormap); plotNets(ROI2Net, NetLB, 10, 1);colorbar;
% title('Correlation Coefficients (r > 0.1) After Regression','FontName','Calibri','FontSize',14);
% set(gca,'FontName','Calibri','FontSize',10);
% set(gcf,'Color','w');

%% Calculate the amount of static connections before and after regression
% Correlation coefficient thresholds
weak_FC   = 0.2;
medium_FC = 0.5;
strong_FC = 0.8;

%%% Static FC before regression
% Count the number of correlations in each category for matrix 1
weak_count_matrix1 =   sum(fc_matrix_1(:) > weak_FC & fc_matrix_1(:) < medium_FC);
medium_count_matrix1 = sum(fc_matrix_1(:) >= medium_FC & fc_matrix_1(:) < strong_FC);
strong_count_matrix1 = sum(fc_matrix_1(:) >= strong_FC);

% Calculate the total number of correlations in matrix 1
total_correlations_matrix1 = numel(fc_matrix_1);

% Calculate the percentage of each category in matrix 1
weak_percentage_matrix1 =   (weak_count_matrix1 / total_correlations_matrix1) * 100;
medium_percentage_matrix1 = (medium_count_matrix1 / total_correlations_matrix1) * 100;
strong_percentage_matrix1 = (strong_count_matrix1 / total_correlations_matrix1) * 100;

%%% Static FC after regression
% Count the number of correlations in each category for matrix 2
weak_count_matrix2 = sum(fc_matrix_2(:) > weak_FC & fc_matrix_2(:) < medium_FC);
medium_count_matrix2 = sum(fc_matrix_2(:) >= medium_FC & fc_matrix_2(:) < strong_FC);
strong_count_matrix2 = sum(fc_matrix_2(:) >= strong_FC);

% Calculate the total number of correlations in matrix 2
total_correlations_matrix2 = numel(fc_matrix_2);

% Calculate the percentage of each category in matrix 2
weak_percentage_matrix2 =   (weak_count_matrix2 / total_correlations_matrix2) * 100;
medium_percentage_matrix2 = (medium_count_matrix2 / total_correlations_matrix2) * 100;
strong_percentage_matrix2 = (strong_count_matrix2 / total_correlations_matrix2) * 100;

%%% Calculate the difference before and after regression
weak_percentage_difference   = weak_percentage_matrix2 - weak_percentage_matrix1;
medium_percentage_difference = medium_percentage_matrix2 - medium_percentage_matrix1;
strong_percentage_difference = strong_percentage_matrix2 - strong_percentage_matrix1;

%%%% Overall difference in connections > 0.2
connections_above_0_2_matrix1 = sum(fc_matrix_1(:) > weak_FC);
perc_above_0_2_matrix1 = (connections_above_0_2_matrix1 / total_correlations_matrix1) * 100;
connections_above_0_2_matrix2 = sum(fc_matrix_2(:) > weak_FC);
perc_above_0_2_matrix2 = (connections_above_0_2_matrix2 / total_correlations_matrix2) * 100;
overall_decrease_connections_0_2 = perc_above_0_2_matrix1 - perc_above_0_2_matrix2;

fprintf('Percentage of weak correlations in matrix 2: %.2f%% (Difference: %.2f%%)\n', weak_percentage_matrix2, weak_percentage_difference);
fprintf('Percentage of medium correlations in matrix 2: %.2f%% (Difference: %.2f%%)\n', medium_percentage_matrix2, medium_percentage_difference);
fprintf('Percentage of strong correlations in matrix 2: %.2f%% (Difference: %.2f%%)\n', strong_percentage_matrix2, strong_percentage_difference);
fprintf('Overall decrease of connections > 0.2 in matrix 2 compared to matrix 1: %.2f%%\n', overall_decrease_connections_0_2);