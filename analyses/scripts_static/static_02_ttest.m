% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESIS TESTING: STATIC FC DIFFERENCES BETWEEN GROUPS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear; clc;
dirhead='...\Analysis';
%% Get the data
% Load the averaged and subject static FC matrices
datadir=[dirhead, '\final_results\Static_FC\'];
fc_matrix_ASD_average = load([datadir, 'averaged_correlation_matrix_ASD.mat']);
fc_matrix_TD_average  = load([datadir, 'averaged_correlation_matrix_TD.mat']);
fc_matrix_ASD = load([datadir, 'subjects_correlation_matrix_ASD.mat']);
fc_matrix_TD  = load([datadir, 'subjects_correlation_matrix_TD.mat']);

% Access the loaded matrix inside the structure
fc_matrix_ASD_average = fc_matrix_ASD_average.averageCorrelationMatrix;
fc_matrix_TD_average = fc_matrix_TD_average.averageCorrelationMatrix;
fc_matrix_ASD = fc_matrix_ASD.correlationMatrices;
fc_matrix_TD = fc_matrix_TD.correlationMatrices;

% Get atlas paramaters as in Xu et al. (2023)
AtlasTable = readtable([dirhead, '\data_static\resources\QPP_atlas.csv']);
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
%% Compute uncorrected p-values
numParcels = size(fc_matrix_TD,1); % TD and ASD matrices have the same size
p_values = zeros(numParcels, numParcels);
for i = 1:numParcels
    for j = 1:numParcels
        % Extract the 21 values from each cell
        values1 = squeeze(fc_matrix_ASD(i, j, :));
        values2 = squeeze(fc_matrix_TD(i, j, :));
        % Perform two-sample t-test
        [~, p] = ttest2(values1, values2);
        % Store the p-value
        p_values(i, j) = p;
    end
end
clear values1;clear values2;clear i;clear j;clear table_ef;

% Set the diagonal values to 1 (this is necessary because variance is 0 for
% identical values in the ttest, outputting NaNs)
for i = 1:numParcels
    p_values(i,i) = 1;
end
%% FDR using the method introduced in Storey (2001)
R = 0.05:0.01:0.95;
[pFDR, Q, aprioriProb] = mafdr(p_values(:),'Lambda',R,'Method','bootstrap','Showplot',true);
clear R;
% pFDR represents the positive FDR. 
% Q represents q-values which are measures of hypothesis testing error for all observations in P-values (input).
% aprioriProb represents the a priori probability that the null hypothesis is true.

% Reshape the pFDR and q-values vector back to the 246x246 matrix
q_values = reshape(Q, 246, 246);
pFDR = reshape(Q, 246, 246);
%% Plotting
%%% Plot p-values
figure;
imagesc(p_values,[0 1]);colorbar;plotNets(ROI2Net,NetLB,10,1);colormap(jet)       
title('P-values Before Regression','FontName','Calibri','FontSize',14);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w');
f = gcf;
plotname1 = 'Pvalues_before_reg.png';
exportgraphics(f,[dirhead, '\final_results\Static_FC\figs\', plotname1],'Resolution',300)
plotname2 = 'Pvalues_before_reg.fig';
savefig(f,[dirhead, '\final_results\Static_FC\figs\', plotname2]);

%%% Plot the histogram of p-values
figure;
histogram(p_values,'FaceColor','blue','FaceColor','blue','LineStyle','none');
box off; grid off;
xlabel('P-value','FontName','Calibri','FontSize',12);
ylabel('Frequency','FontName','Calibri','FontSize',12);
title('Histogram of P-values Before Regression','FontName','Calibri','FontSize',14);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w');
f = gcf;
plotname1 = 'Pvalues_hist_before_reg.png';
exportgraphics(f,[dirhead, '\final_results\Static_FC\figs\', plotname1],'Resolution',300)
plotname2 = 'Pvalues_hist_before_reg.fig';
savefig(f,[dirhead, '\final_results\Static_FC\figs\', plotname2]);

%%% Plot q-values
figure;
imagesc(q_values,[0 1]);colorbar;plotNets(ROI2Net,NetLB,10,1);colormap(jet)       
title('Q-values Before Regression','FontName','Calibri','FontSize',14);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w');
f = gcf;
plotname1 = 'Qvalues_before_reg.png';
exportgraphics(f,[dirhead, '\final_results\Static_FC\figs\', plotname1],'Resolution',300)
plotname2 = 'Qvalues_before_reg.fig';
savefig(f,[dirhead, '\final_results\Static_FC\figs\', plotname2]);

%%% Plot q values that survive FDR < 0.05
%%%% Define significance threshold & values below it
significance_level = 0.05;
q_values_thresholded = q_values;
q_values_thresholded(q_values >= significance_level) = NaN;
%%%% Plot only the values below the threshold
figure;
imagesc(q_values_thresholded,[0 significance_level]); % Only FDR < 0.05
custom_colormap = [0.5 0.5 0.5; jet(256)]; % Add gray to the colormap
colormap(custom_colormap); plotNets(ROI2Net, NetLB, 10, 1);colorbar;
title('Q-values (FDR < 0.05) Before Regression','FontName','Calibri','FontSize',14);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w');
f = gcf;
plotname1 = 'Qvalues_sig_before_reg.png';
exportgraphics(f,[dirhead, '\final_results\Static_FC\figs\', plotname1],'Resolution',300)
plotname2 = 'Qvalues_sig_before_reg.fig';
savefig(f,[dirhead, '\final_results\Static_FC\figs\', plotname2]);
%% Save output
save([dirhead, '\final_results\Static_FC\Pvalues_matrix.mat'],"p_values");
save([dirhead, '\final_results\Static_FC\Qvalues_matrix.mat'],"q_values");
save([dirhead, '\final_results\Static_FC\Qvalues_sig_matrix.mat'],"q_values_thresholded");
save([dirhead, '\final_results\Static_FC\pFDR_matrix.mat'],"pFDR");