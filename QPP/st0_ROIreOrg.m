%% Setup parameters
%%% load unreorganized data
dataIN='./Input/preprocessed_combined_TD.mat'; D=load(dataIN); D0=D.combined; 
%%% define output filename
dataOUT='./Input/GU_TD_SUB_gsr.mat'; 
%%% load parcel label and network spreedsheet
AtlasTable = readtable(['./resources/QPP_atlas.csv']);
Label='Label';         % this should be variable name of the numerical label of ROIs
System='network7_shortname'; % this is the variable name of the numerical label of functional networks
%%% add function files
addpath('./QPPfv0922/')
%% code
AtlasTable.Label=eval(['AtlasTable.' Label]); AtlasTable.System=eval(['AtlasTable.' System]);
AtlasTable = sortrows(AtlasTable,'Label','ascend');
[NetLB, net_index, ROI2Net]=unique(AtlasTable.System,'stable'); 
[nsbj, nscn]=size(D0); [nroi,ntimepoint]=size(D0{1,1});
if ~isfield(D, 'MotionInf')
    for i=1:nsbj, for j=1:nscn,[~,ntimepoint]=size(D0{i,j});  MotionInf{i,j}=[1:ntimepoint]; end; end
end
save(dataOUT,'D0','MotionInf','ROI2Net','NetLB');