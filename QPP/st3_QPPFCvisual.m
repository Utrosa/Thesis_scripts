clear; clc; close all
%% Setup data input & visualization parameters
% dataext='HumanVisual_allvols'; % extended filename=[data '_' ext];
dataext='GU_TD_SUB_gsr'; % extended filename=[data '_' ext];
p2param=['Params_' dataext '.mat']; load(p2param);

runM=2; % QPP running method
% runM: 1 -concatenate all D{i,j} as a whole group and detect group QPP
%       2 -concatenate all D{i,:} and detect QPP from all scans of each subject
%       3 -concatenate all D{:,j} and detect QPP from all subjects of each
%       scan
rbstScrn=1; % control for robust QPP detection
%rbstScrn:  1 - scan all possible initial segments for robust QPP detection
%           0 - scan randomly slected initial segments for fast QPP detection
Pselect=[1];  %select QPP#s to be visualized
Gselect=[1:4];  %select groups to be compared
bin=0.08;
%% Automatically load data & other hidden parameters
fprintf('QPP & FC Visualizations...\n'); 
addpath(p2qppf);
d2O='./Output/';  % directory to outputs files
if runM==1
    Ng=1;Gselect=1; a0=[d2O 'GrpQPP/']; if ~exist(a0,'dir'), mkdir(a0); end % dir2 save GrpQPPs  
    a=[a0 'interm/']; if ~exist(a,'dir'), mkdir(a); end % save intermediate files
    indn='Grp'; 
elseif runM==2
    Ng=length(Gselect); a0=[d2O 'SbjQPP/']; if ~exist(a0,'dir'), mkdir(a0); end % dir2 save SbjQPPs  
    a=[a0 'interm/']; if ~exist(a,'dir'), mkdir(a); end % save intermediate files
    indn='Sbj';
elseif runM==3   
    Ng=length(Gselect); a0=[d2O 'ScanQPP/']; if ~exist(a0,'dir'), mkdir(a0); end % dir2 save ScanQPPs  
    a=[a0 'interm/']; if ~exist(a,'dir'), mkdir(a); end % save intermediate files
    indn='Scn';
else, error('Unidenfied value for runM.\n')
end
p2S=cell(Ng,nP); p2S0=cell(Ng,1); % pth2 save SbjQPPs
for ig=Gselect % gp id
    for ip=Pselect  % QPP id
        p2S{ig,ip}=[a dataext '_' indn num2str(ig) '_rbst' num2str(rbstScrn) '_qpp' num2str(ip)];
    end
    p2S0{ig}=[a0 dataext '_' indn num2str(ig) '_rbst' num2str(rbstScrn) '_QPPs'];
     
end
%% Visualization
ct=0; nPs=length(Pselect);
for ip=Pselect 
    for ig=Gselect
        load(p2S0{ig},'QPPs','TMXs','Cs','METs', 'QPPas', 'TMXas', 'METas', 'Cas', ...
        'Ds', 'Drs', 'Crs', 'FCrs', 'ROI2Net', 'iROI2Net','NetLB', ...
        'ntlist','ITP','PL','ssg','tres',...        
        'cth13','cth45','paramQPPf2','paramQPPf4');      
 
        ct=ct+1;
        f1=figure(1); %qpps
        subplot(nPs,Ng,ct); imagesc(QPPs{ip,1}(iROI2Net,:),[-1 1]); colormap(jet)       
        plotNets(ROI2Net,NetLB, PL(ip),0); colorbar
        if ig==1, ylabel(['QPP #' num2str(ip)]); end
        if ip==1, title([indn num2str(ig)],'FontSize', 8); end

        f2=figure(2); %reverse qpps
        subplot(nPs,Ng,ct); imagesc(QPPs{ip,2}(iROI2Net,:),[-1, 1]); colormap(jet)
        plotNets(ROI2Net,NetLB, PL(ip),0); colorbar
        if ig==1, ylabel(['rph QPP #' num2str(ip)]); end
        if ip==1, title(['reverse phase (rph) QPP: ' indn num2str(ig)],'FontSize', 8); end
              
        f3=figure(3); %qpp sliding correlations        
        subplot(nPs,Ng,ct); C1=Cs(ip,:); 
        TMX1=TMXs{ip,1}; Met1=METs{ip,1}; TMX2=TMXs{ip,2}; Met2=METs{ip,2};
        plot(C1,'b'); hold on; plot(TMX1,C1(TMX1),'bv'); plot(TMX2,C1(TMX2),'bo'); 
        axis([0 length(C1) -1 1]); 
        set(gca,'XTick',1:sum(ntlist)/4:sum(ntlist),'YTick',-1:0.2:1); grid on
        title({['sliding correlation: ' indn num2str(ig)],['median max: ' num2str(0.01*round(100*Met1(1))) ...
        ', median \Deltat_{max}:' num2str(0.1*round(10*Met1(2))) ' tps' ...
        ', #max:' num2str(Met1(3))], ['median min: ' num2str(0.01*round(100*Met2(1))) ...
        ', median \Deltat_{min}:' num2str(0.1*round(10*Met2(2))) ' tps' ...
        ', #min:' num2str(Met2(3))]},'FontSize', 8,'fontweight','normal');
        if ig==1,  ylabel(['QPP #'  num2str(ip)],'FontSize', 8); end
           
        f5=figure(4); %qpp sliding correlations before and after qpp regression     
        subplot(nPs,Ng,ct); Cr2=Crs{ip}(ip,:);
        plot(C1,'b'); hold on; plot(Cr2,'r'); axis([0 length(C1) -1 1]); 
        set(gca,'XTick',1:sum(ntlist)/4:sum(ntlist),'YTick',-1:0.2:1); grid on          
        legend('before qpp regression','after ~');
        if ig==1,  ylabel(['QPP #'  num2str(ip)],'FontSize', 8); end
        if ip==1,  title({'sliding corr', ['timecourse: ' indn num2str(ig)]},'FontSize', 8); end
        
        f6=figure(5); %histograms of sliding correlation timecourse before & after qpp regression     
        subplot(nPs,Ng,ct); 
        hist2distns(C1, Cr2,'before qpp reg', 'after ~', 'b','r', bin);
        if ig==1,  ylabel(['QPP #'  num2str(ip)],'FontSize', 8); end
        legend('before qpp regression','after ~'); 
        if ip==1, title(['sliding corr histogram: ' indn num2str(ig)],'FontSize', 8); end
    
        f7=figure(6); %FC matrix before and after qpp regression     
        subplot(nPs,Ng,ct); FC=corr(Ds{ip}(iROI2Net,:)'); FC1=triu(FC,0); FC2=tril(FCrs{ip},-1);
        imagesc(FC1+FC2,[-1 1]); colorbar; plotNets(ROI2Net,NetLB, PL(ip),1); colormap(jet)       
        if ig==1,  ylabel(['QPP #'  num2str(ip)],'FontSize', 8); end
        if ip==1, title({['FC before (up) & after (low)'], ['qpp reg: ' indn num2str(ig)]},'FontSize', 8); end
    
    end   
end