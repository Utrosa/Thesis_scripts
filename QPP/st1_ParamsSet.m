clear; clc;  close all
%% For saving & reading path-names
data='GU_TD_SUB_gsr'; % the filename of the data file
p2data=['./Input/' data '.mat']; % path of data file which must includes following parameters:
% D0:   a nsbj X nscn cell matrix. Each cell has a nroi X ntimepoints 
%       matrix of EPI timeseries
% MotionInf: a nsbj X nscn cell matrix. Each cell >=1 segments of
%       timepoints without significant motions.
% ROI2Net:  a (nroi X 1) vector including the network index of each ROI
% NetLB: a (nnet X 1) cell vector including the label for each network
% TR: repetition time for each subject

ext=['']; % filename extension for the parameter file
if ~isempty(ext), p2param=['Params_' data '_' ext]; else, p2param=['Params_' data];  end% parameter filename
p2qppf='./QPPfv0922/'; %directory of QPP functions
%% Set up global parameters 
%%%%%% for QPP detection %%%%%%
% PL: a (nP X 1) vector of QPP widow length, ~20s in humans. E.g.,
% PL(ip)=20s/TR.
PL=[10];
%% For QPP detection: QPPf1detectRbs.m
% cth13 & cth45: correlation threshold for QPP1-QPP3 & for QPP4-QPP5;
cth13=[0.1, 0.2]; cth45=[0.1, 0.2];

%% For QPP phase adjustment: QPPf2phadj.m
% cthph: similarity threshold when phase-adjusting (phadj) a QPP
% s: 1 or 0
%    1 -strict phase adjustment
%    0 -relaxed phase adjustment
cthph=0.88; s=0; 

% Reference parcel(s) ID for QPP phase adjustment: 
% One can >=1 parcel(s). These parcels will have QPP waveform starting from 
% rising postive values (e.g. sine wave).
nP=length(PL); sdph=cell(nP,1); 
sdph{1}={163;18}; sdph{2}={163; 87}; sdph{3}={18;87}; sdph(4:5)={{18}}; 
%% For network analysis (e.g., QPPf4regscnRbst.m)
% fz: 1 or 0. 
%     1 -output matrix FCr will be the pearson correlation matrix
%     0 -output matrix FCr will be the Fisher Z-Transformation of the pearson correlation
fz=0; 
save([p2param '.mat']);