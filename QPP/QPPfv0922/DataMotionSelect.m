function [Ddm, ntlist]=DataMotionSelect(D, MotionInf)
%% Inputs & Outputs
%%Input:
% D: a nsbj X nscn cell matrix of original EPI timeseries.
% MotionInf: a nsbj X nscn cell matrix, each cell includes the non-moving timepoints

%%Output:
% Ddm: a nROI X (nsbj*nscn) matrix of cancatenated EPI timeseries without
% significant motions
% ntlist: a list of timepoint number for each scan
%% 
Ddm=[];  ntlist=[]; scan_ct=0;
[nsbj,nscn]=size(MotionInf);
for isbj=1:nsbj
    for iscn=1:nscn        
        Dtmp=D{isbj,iscn}; timeSelect=MotionInf{isbj,iscn};
        if isempty(timeSelect)
            continue;
        elseif size(timeSelect,2)>2            
            scan_ct=scan_ct+1;
            Dtmp=Dtmp(:,timeSelect); nt=size(Dtmp,2); 
            Dtmp=zscore(Dtmp,[],2);
            Ddm=cat(2,Ddm,Dtmp); ntlist=[ntlist;nt];      

        elseif size(timeSelect,2)==2 % two segments
            for seg=1:size(timeSelect,2)                    
                tmp=timeSelect{seg}; 
                Dtmp1=Dtmp(:,tmp); scan_ct=scan_ct+1;
                nt=size(Dtmp1,2); ntlist=[ntlist;nt];
                Dtmp1=zscore(Dtmp1,[],2);
                Ddm=cat(2,Ddm,Dtmp1); 
            end
        end        
    end
end

    