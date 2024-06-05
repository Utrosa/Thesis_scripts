% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% STATIC FUNCTIONAL CONNECTIVITY AFTER REGRESSION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear; clc; close all;
%% Define type of analysis: scan QPP or subject QPP
cs_type = 'scan';
%cs_type = 'subject'; % Not recommended
%% Get atlas parameters as in Xu et al. (2023)
dirhead = '...\Analysis';
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
clear System; clear Label; clear a;
%% Load the Regressed Functional Connectivity Matrices
%%% Load structs with FCrs QPP output. ASD_files & TD_files are created
%%% with the QPP_02_group_level_comparison.m script.
data_ASD = load([dirhead,'\final_results\Dynamic_FC\ASD_files.mat']);
data_TD = load([dirhead,'\final_results\Dynamic_FC\TD_files.mat']);

%%% Extract FCrs data matrices
if cs_type == "subject"
    ASD_FCrs = zeros(246,246,90); % 90 participants across collections
    TD_FCrs  = zeros(246,246,90);
    for i = 1:length(data_TD.TD_files) % The same length for ASD and TD
        TD_FCrs(:,:,i) = data_TD.TD_files(i).data.FCrs{1,1};
        ASD_FCrs(:,:,i) = data_ASD.ASD_files(i).data.FCrs{1,1};
    end
    % Average regressed FC matrices across subjects for each group
    FC_REG_ASD = mean(ASD_FCrs, 90);
    FC_REG_TD  = mean(TD_FCrs, 90);

elseif cs_type == "scan"
    ASD_FCrs = zeros(246,246,3); % 3 data collections
    TD_FCrs  = zeros(246,246,3);
    for i = 1:length(data_TD.TD_files) % The same length for ASD and TD
        TD_FCrs(:,:,i) = data_TD.TD_files(i).data.FCrs{1,1};
        ASD_FCrs(:,:,i) = data_ASD.ASD_files(i).data.FCrs{1,1};
    end
    % Average regressed FC matrices across data collections
    FC_REG_ASD = mean(ASD_FCrs, 3);
    FC_REG_TD  = mean(TD_FCrs, 3);
end

%%% Save the regressed FC matrices
save([dirhead,'\final_results\Static_FC\subjects_correlation_matrix_reg_ASD.mat'],"ASD_FCrs");
save([dirhead,'\final_results\Static_FC\subjects_correlation_matrix_reg_TD.mat'], "TD_FCrs");
save([dirhead,'\final_results\Static_FC\averaged_correlation_matrix_reg_ASD.mat'],"FC_REG_ASD");
save([dirhead,'\final_results\Static_FC\averaged_correlation_matrix_reg_TD.mat'], "FC_REG_TD");
%% Plotting
%%%% Averaged Static FC Matrix for Control participants
figure;
imagesc(FC_REG_TD,[-1 1]);colorbar;plotNets(ROI2Net,NetLB,10,1);colormap(jet)   
title('Static FC After Regression: TD','FontName','Calibri','FontSize',14);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w');
f = gcf;
exportgraphics(f,[dirhead, '\final_results\Static_FC\figs\staticFC_after_reg_TD.png'],'Resolution',300)
savefig(f,[dirhead, '\final_results\Static_FC\figs\staticFC_after_reg_TD.fig']);

%%%% Averaged Static FC Matrix for Autistic participants
figure;
imagesc(FC_REG_ASD,[-1 1]);colorbar;plotNets(ROI2Net,NetLB,10,1);colormap(jet)   
title('Static FC After Regression: ASD','FontName','Calibri','FontSize',14);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w');
f = gcf;
exportgraphics(f,[dirhead, '\final_results\Static_FC\figs\staticFC_after_reg_ASD.png'],'Resolution',300)
savefig(f,[dirhead, '\final_results\Static_FC\figs\staticFC_after_reg_ASD.fig']);
%% Perform multiple two-sample t tests
fc_matrix_TD  = TD_FCrs;
fc_matrix_ASD = ASD_FCrs;
% Compute uncorrected p-values
numParcels = size(fc_matrix_TD,1); % ASD and TD matrices have the same size
p_values = zeros(numParcels, numParcels);
for i = 1:numParcels
    for j = 1:numParcels
        % Extract the values across the third dimension
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
% pFDR represents the positive FDR. 
% Q represents q-values which are measures of hypothesis testing error for
% all observations in P-values (input).
% aprioriProb represents the a priori probability that the null hypothesis
% is true.
R = 0.05:0.01:0.95;
[pFDR, Q, aprioriProb] = mafdr(p_values(:),'Lambda',R,'Method','bootstrap','Showplot',true);
clear R;

% Reshape the pFDR and q-values vector back to the 246x246 matrix
q_values = reshape(Q, 246, 246);
pFDR = reshape(Q, 246, 246);
%% Plotting
%%% Plot p-values
figure;
imagesc(p_values,[0 1]);colorbar;plotNets(ROI2Net,NetLB,10,1);colormap(jet)       
title('P-values After Regression','FontName','Calibri','FontSize',14);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w');
f = gcf;
plotname1 = 'Pvalues_after_reg.png';
exportgraphics(f,[dirhead, '\final_results\Static_FC\figs\', plotname1],'Resolution',300)
plotname2 = 'Pvalues_after_reg.fig';
savefig(f,[dirhead, '\final_results\Static_FC\figs\', plotname2]);

%%% Plot the histogram of p-values
figure;
histogram(p_values,'FaceColor','blue','EdgeColor','black'); box off;
xlabel('P-value','FontName','Calibri','FontSize',12);
ylabel('Frequency','FontName','Calibri','FontSize',12);
title('Histogram of P-values After Regression','FontName','Calibri','FontSize',14);
set(gca,'FontName','Calibri','FontSize',10);grid off;
set(gcf,'Color','w');
f = gcf;
plotname1 = 'Pvalues_hist_after_reg.png';
exportgraphics(f,[dirhead, '\final_results\Static_FC\figs\', plotname1],'Resolution',300)
plotname2 = 'Pvalues_hist_after_reg.fig';
savefig(f,[dirhead, '\final_results\Static_FC\figs\', plotname2]);

%%% Plot q-values
figure;
imagesc(q_values,[0 1]);colorbar;plotNets(ROI2Net,NetLB,10,1);colormap(jet)       
title('Q-values After Regression','FontName','Calibri','FontSize',14);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w');
f = gcf;
plotname1 = 'Qvalues_after_reg.png';
exportgraphics(f,[dirhead, '\final_results\Static_FC\figs\', plotname1],'Resolution',300)
plotname2 = 'Qvalues_after_reg.fig';
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
title('Q-values (FDR < 0.05) After Regression','FontName','Calibri','FontSize',14);
set(gca,'FontName','Calibri','FontSize',10);
set(gcf,'Color','w');
f = gcf;
plotname1 = 'Qvalues_sig_after_reg.png';
exportgraphics(f,[dirhead, '\final_results\Static_FC\figs\', plotname1],'Resolution',300)
plotname2 = 'Qvalues_sig_after_reg.fig';
savefig(f,[dirhead, '\final_results\Static_FC\figs\', plotname2]);
%% Save output
save([dirhead, '\final_results\Static_FC\Pvalues_reg_matrix.mat'],"p_values");
save([dirhead, '\final_results\Static_FC\Qvalues_reg_matrix.mat'],"q_values");
save([dirhead, '\final_results\Static_FC\Qvalues_reg_sig_matrix.mat'],"q_values_thresholded");
save([dirhead, '\final_results\Static_FC\pFDR_reg_matrix.mat'],"pFDR");