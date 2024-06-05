% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% STATIC FUNCTIONAL CONNECTIVITY MATRICES: Correlation Coefficients
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear;clc;close all;
% Set the parameters
dirhead = '...\Analysis';
%group = 'ASD';
group = 'TD';
%% Get all the prerequisites
% Get a list of subjects: sublist_groupname.txt
subs_path = sprintf('%s\\data_static\\sublist_%s.txt', dirhead, group);
subs_list = fopen(subs_path, 'r');
subs = fscanf(subs_list, "%c");
fclose(subs_list);
subjects = regexp(subs, 'sub-\d{5}', 'match');
clear subs; clear subs_list; clear subs_path;

% Get the data
datadir = [dirhead, '\data_static\Input\'];
alldata = cell(numel(subjects),1);
for i = 1:numel(subjects)
    subjectdir = [datadir, subjects{i}];
    filepath   = fullfile(subjectdir, sprintf('%s_wmcsfr_gsr_fil_on_parcel_fan246.mat', subjects{i}));
    loadedfile = load(filepath);
    ts_parcel  = loadedfile.ts_zscore; % Load z-scored data
    ts_parcel(1,:) = [];
    alldata{i} = ts_parcel;
end
clear ts_parcel; clear loadedfile; clear filepath; clear subjectdir;
clear datadir;

% Get atlas parameters as in Xu et al. (2023)
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
clear i; clear Label; clear System; clear loadedfile;
%% Correlation matrices
correlationMatrices = zeros(246,246,numel(alldata));
%correlationMatrices = struct();
for i = 1:numel(alldata)

    % Get EPI timeseries for each subject & number of brain regions
    meanTimeCourse = alldata{i};
    numParcels = size(meanTimeCourse, 1);

    % Calculate Pearson correlation coefficients between each pair of ROIs
    reorderedTs = meanTimeCourse(iROI2Net,:);
    
    % Transpose: rows are observations (time), columns variables (ROIs)
    correlationMatrix = corrcoef(reorderedTs');

    % Fisher's z-transform
    z_scored_correlationMatrix = atanh(correlationMatrix);

    % Correct the diagonal values. When atanh() is applied to a value of 1 
    % the result is mathematically undefined, so MATLAB returns Inf. The 
    % diagnoal values: r = 1.
    z_scored_correlationMatrix(1:numParcels+1:end) = 1;

    % Save a matrix for each subject
    correlationMatrices(:,:,i) = z_scored_correlationMatrix;
    % fieldName = sprintf('subject_%d', i);
    % correlationMatrices.(fieldName) = correlationMatrix;

end 
output1 = sprintf('subjects_correlation_matrix_%s.mat', group);
save([dirhead, '\final_results\Static_FC\', output1], 'correlationMatrices');

% Average correlation matrices across subjects
averageCorrelationMatrix = mean(correlationMatrices, 3);
output2 = sprintf('averaged_correlation_matrix_%s.mat', group);
save([dirhead, '\final_results\Static_FC\', output2], 'averageCorrelationMatrix');
%% Plot the Averaged Static FC Matrix: Correlation (per group)
figure;
imagesc(averageCorrelationMatrix,[-1 1]); colorbar; 
plotNets(ROI2Net,NetLB,10,1); colormap(jet)   
corr_title = sprintf('Static FC Before Regression: %s', group);
title(corr_title,'FontName', 'Calibri','FontSize',14);
set(gca, 'FontName', 'Calibri', 'FontSize', 10);
set(gcf, 'Color', 'w');
f = gcf;
output3 = sprintf('staticFC_before_reg_%s.png', group);
exportgraphics(f,[dirhead, '\final_results\Static_FC\figs\', output3],'Resolution',300)

output4 = sprintf('staticFC_before_reg_%s.fig', group);
savefig(f,[dirhead, '\final_results\Static_FC\figs\', output4]);