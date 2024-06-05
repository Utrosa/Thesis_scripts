clear; clc; close all
%% Setup data input, output directory & running method
% dataext='HumanVisual_test'; % extended filename=[data '_' ext];
dataext='GU_TD_SUB_gsr'; % extended filename=[data '_' ext];
runM=1; % QPP running method
% runM: 1 -concatenate all D{i,j} as a whole group and detect group QPP
%       2 -concatenate all D{i,:} and detect QPP from all scans of each subject
%       3 -concatenate all D{:,j} and detect QPP from all subjects of each
%       scan
rbstScrn=1; % control for robust QPP detection
%rbstScrn:  1 - scan all possible initial segments for robust QPP detection
%           0 - scan randomly slected initial segments for fast QPP detection
%% Automatically load data & other hidden parameters
fprintf('Loading data & hidden parameters\n'); 
p2param=['Params_' dataext '.mat']; load(p2param); addpath(p2qppf);
load(p2data, 'D0','MotionInf','ROI2Net','NetLB'); [nsbj,nscn]=size(D0); 

ssg=ones(nP,1); ssg(2:end)=PL(2:end); tres=0.7; 
paramQPPf2=param_QPPf2(nP, PL, cthph, sdph, s);
% paramQPPf4=param_QPPf4(nP, PL, ibY, iROI2Net, fz); 
paramQPPf4=param_QPPf4(nP, PL, ROI2Net, fz);iROI2Net= paramQPPf4.iROI2Net; 
ITPstp=20; % step to show progress of algorithm when running ITP times

d2O='./Output/'; if ~exist(d2O,'dir'), mkdir(d2O); end % directory to outputs files
if runM==1
    Ng=1; a0=[d2O 'GrpQPP/']; if ~exist(a0,'dir'), mkdir(a0); end % dir2 save GrpQPPs  
    a=[a0 'interm/']; if ~exist(a,'dir'), mkdir(a); end % save intermediate files
    indn='Grp'; 
elseif runM==2
    Ng=nsbj; a0=[d2O 'SbjQPP/']; if ~exist(a0,'dir'), mkdir(a0); end % dir2 save SbjQPPs  
    a=[a0 'interm/']; if ~exist(a,'dir'), mkdir(a); end % save intermediate files
    indn='Sbj';
elseif runM==3   
    Ng=nscn; a0=[d2O 'ScanQPP/']; if ~exist(a0,'dir'), mkdir(a0); end % dir2 save ScanQPPs  
    a=[a0 'interm/']; if ~exist(a,'dir'), mkdir(a); end % save intermediate files
    indn='Scn';
else, error('Unidenfied value for runM.\n')
end
p2S=cell(Ng,nP); p2S0=cell(Ng,1); % pth2 save SbjQPPs
for ig=1:Ng % gp id
    for ip=1:nP % QPP id
        p2S{ig,ip}=[a dataext '_' indn num2str(ig) '_rbst' num2str(rbstScrn) '_qpp' num2str(ip)];
    end
    p2S0{ig}=[a0 dataext '_' indn num2str(ig) '_rbst' num2str(rbstScrn) '_QPPs'];
     
end
%% Computation
delete(gcp('nocreate'));
myCluster = parcluster('local'); myCluster.NumWorkers = 56;  % 'Modified' property now TRUE
saveProfile(myCluster); 
parpool(10,'IdleTimeout',100000000)
for ig=1:Ng
    if runM==1,     D00=D0; MotionInf1=MotionInf;
    elseif runM==2, D00=D0(ig,:); MotionInf1=MotionInf(ig,:);
    elseif runM==3, D00=D0(:,ig); MotionInf1=MotionInf(:,ig); 
    end    
     
    [D, ntlist]=DataMotionSelect(D00, MotionInf1); nscng=length(ntlist); D1=D;
    paramQPPf1=param_QPPf1(nP, ntlist, PL, cth13, cth45, ssg);
    QPPs=cell(nP,2);  TMXs=QPPs; METs=QPPs;%% QPP templates and their reverse templates    
    QPPas=cell(nP,1); TMXas=QPPas; METas=QPPas; Cas=QPPas; FCrs=QPPas; 
    Cs=zeros(nP,sum(ntlist),'single'); % sliding window correlation    
    Ds=cell(nP,1); Drs=Ds; Crs=Ds; 
     
    for ip=1:nP       
        tic; if ip~=1, D1=load(p2S{ig,ip-1},'Dr'); D1=D1.Dr; end
        
        isip=sprintf([indn '%d-QPP%d-'], ig,ip); 
        if rbstScrn==1
            ITP=paramQPPf1.ITProbust{ip}; % for robust detection
        else
            ITP=paramQPPf1.ITPfast{ip}; % for fast detection
        end
        
        [QPP,TMX,C,MET,ITER,TMPL,TMXTMPL,CTMPL,SCMX]=QPPf1detectRbst(D1,nscng,ntlist, ...
                        paramQPPf1.PL(ip), ...
                        paramQPPf1.cth{ip}, ...
                        paramQPPf1.ncth1(ip), ...
                        paramQPPf1.nitr, ...
                        paramQPPf1.ssg(ip),...
                        ITP, ...
                        paramQPPf1.PLh{ip},...
                        tres, [isip 'f1detect'],ITPstp);
                                    
        fprintf([isip 'f2phadj\n']);
        [QPPa,TMXa,Ca,METa,SDa,SD,flga,cT1Tj,nsim]=QPPf2phadj(QPP,TMPL,SCMX,TMXTMPL,CTMPL, ...
                        paramQPPf2.tsh(ip), ...
                        paramQPPf2.PLc{ip}, ...
                        paramQPPf2.cthph, ...
                        paramQPPf2.sdph{ip}, ...
                        paramQPPf2.s,tres);

        fprintf([isip 'f3xresid\n']);
        if ip~=1
            QPP0=QPP; QPPa0=QPPa; C0=C; SD0=SD; SD=nan; SDa0=SDa; SDa=nan;
            [QPP,C,MET(1),QPPa,Ca,METa(:,1)]=QPPf3xresidRbst...
                (D,ntlist, TMX,TMXa,paramQPPf1.PLh{ip},paramQPPf2.PLc{ip},nscng,ssg(ip));
        end
        
        fprintf([isip 'f3detect_reverse\n']);
        [QPPn, TMXn, METn]=QPPf1detectN(D, C, ntlist, PL(ip), paramQPPf1.cth{ip},tres);

        QPPs{ip,1}=QPP; QPPs{ip,2}=QPPn; Cs(ip,:)=C;
        fprintf([isip 'f4regscn\n']);
        for ip2=1:ip-1
            if cellfun(@isempty,(QPPs(ip2,1)))
                fl=load(p2S{ig,ip2},'QPP','QPPn','C','TMX','MET','TMXn','METn','QPPa','Ca','TMXa','METa','D','Dr','Cr','FCr'); 
                QPPs{ip2,1}=fl.QPP; QPPs{ip2,2}=fl.QPPn; 
                TMXs{ip2,1}=fl.TMX; TMXs{ip2,2}=fl.TMXn; 
                METs{ip2,1}=fl.MET; METs{ip2,2}=fl.METn; 
                Cs(ip2,:)=fl.C; QPPas{ip2,1}=fl.QPPa; 
                TMXas{ip2,1}=fl.TMXa; METas{ip2,1}=fl.METa; Cas{ip2}=fl.Ca;
                Ds{ip2,1}=fl.D; Drs{ip2,1}=fl.Dr; Crs{ip2}=fl.Cr; FCrs{ip2,1}=fl.FCr;
            end
        end
        [Dr,Cr,FCr]=QPPf4regscnRbst(D,ntlist, QPPs(1:ip,1),Cs(1:ip,:),nscng,PL(ip), ...
                                    paramQPPf4.PLc{ip},...
                                    paramQPPf4.iNetL, ...
                                    paramQPPf4.iROI2NetC, ...
                                    paramQPPf4.fz);                                                               

        save(p2S{ig,ip},'D','QPP','TMX','C','MET','ITER','TMPL','TMXTMPL',...
            'CTMPL','SCMX','QPPn', 'TMXn', 'METn', ...
            'QPPa','TMXa','Ca','METa','SDa','SD','flga',...
            'cT1Tj','nsim','Dr','Cr','FCr', ...
            'ntlist','ITP','PL','ssg','tres',...
            'paramQPPf1','paramQPPf2','paramQPPf4','-v7.3'); 
        
        
        TMXs{ip,1}=TMX; TMXs{ip,2}=TMXn; METs{ip,1}=MET; METs{ip,2}=METn; 
        QPPas{ip,1}=QPPa; TMXas{ip,1}=TMXa; METas{ip,1}=METa; Cas{ip}=Ca;
        Ds{ip,1}=D; Drs{ip,1}=Dr; Crs{ip}=Cr; FCrs{ip,1}=FCr;
        
        if ip~=1, save(p2S{ig,ip},'QPP0','QPPa0','SD0','SDa0','C0','-append'); end
        fprintf([indn '%d-QPP%d %dsec\n\n'],ig,ip,round(toc));

    end
    save(p2S0{ig},'QPPs','TMXs','Cs','METs', 'QPPas', 'TMXas', 'METas', 'Cas', ...
        'Ds', 'Drs', 'Crs', 'FCrs', 'ROI2Net','iROI2Net','NetLB', ...
        'ntlist','ITP','PL','ssg','tres',...        
        'cth13','cth45','paramQPPf2','paramQPPf4','-v7.3');     
end