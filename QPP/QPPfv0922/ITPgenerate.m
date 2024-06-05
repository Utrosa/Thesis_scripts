function [ITP, ITP1]=ITPgenerate(ntlist, PL, ssg, t)
%% Inputs & Outputs
%%Input:
% ntlist: a list of timepoint number for each scan
% PL: the window length of QPP
% ssg: a number for staring segment/scan: ssg=1 for QPP1 detection; ssg=PL
% for QPP2-QPP5 detections

%%Output:
% ITP: all possible initial-segments/sbj for robust QPPdetection
% ITP: randomly selected #initial-segments/scan for fast QPP detection
%% 
ITP=[]; 
for i=1:length(ntlist)
   nt=ntlist(i); nT=sum(ntlist(1:i));
   esg=nt-PL+1; a=nT-nt+(ssg:esg);   
   ITP=[ITP; single(a')];   % for robust detection
end

ITP1=[];
for i=1:length(ntlist)
   nt=ntlist(i);  nT=sum(ntlist(1:i));
   esg=nt-PL+1;  a=nT-nt+(ssg:esg);
   nt2PL=round(nt/PL); nITPps=round(nt2PL./t); % #initial-segments/scan    
   a=a(randperm(length(a))); a=a(1:nITPps); 
   ITP1=[ITP1; single(a')];     
end