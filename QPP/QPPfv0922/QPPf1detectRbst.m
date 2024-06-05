function [QPP,TMX,C,MET,ITER,TMPL,TMXTMPL,CTMPL,SCMX]=...
    QPPf1detectRbst(D,nscn,ntlist, PL,cth,ncth1,nitr,ssg,ITP,PLh,tres,s,stp)
% developed by Nan Xu on June 23, 2021.
[nX,nT]=size(D); nXL=nX*PL; PLe=PL+sum(PLh);

%% Flattening all segments & zscoring
SGf=zeros(nXL,nT,'single'); SGfn=zeros(nXL,nT,'single');
for iscn=1:nscn
    nt=ntlist(iscn); esg=nt-PL+1; 
    nTc=sum(ntlist(1:iscn-1));
for isg=ssg:esg
    t=(isg:isg+PL-1)+nTc; 
    S=D(:,t); S=S(:); SGf(:,nTc+isg)=S;
    S=S-sum(S)/nXL; S=S/sqrt(S'*S); SGfn(:,nTc+isg)=S; 
end
end; clear S

%% Running the QPP algorithm for initial segments of ITP
nITP=length(ITP); TMPL=cell(nITP,1); TMXTMPL=cell(nITP,1); 
CTMPL=zeros(nITP,nT,'single'); SCMX=zeros(nITP,1,'single'); 
ITR=zeros(nITP,1,'single'); e=0.9999;
parfor itp=1:nITP
    c=SGfn(:,ITP(itp))'*SGfn;
    [~,tmx]=findpeaks(c,'MinPeakHeight',cth(1),'MinPeakDistance',PL);
    for iscn=1:nscn
        nt=ntlist(iscn); esg=nt-PL+1; nTc=sum(ntlist(1:iscn-1));
        tmx( tmx==ssg+nTc | tmx==esg+nTc)=[]; 
    end
    nmx=length(tmx);
    cn=c-sum(c)/nT; cn=cn'/sqrt(cn*cn');
    
    cn1=cn; cn2=cn; cn3=cn; itr=1;
    while itr<=nitr 
        if nmx<=1, break; end
        
        T=sum(SGf(:,tmx),2)/nmx; T=T-sum(T)/nXL; T=T/sqrt(T'*T); 
        c=T'*SGfn;
        if itr<=ncth1-1, ith=1; else, ith=2; end   
        [~,tmx]=findpeaks(c,'MinPeakHeight',cth(ith),'MinPeakDistance',PL);
        for iscn=1:nscn
            nt=ntlist(iscn); esg=nt-PL+1; nTc=sum(ntlist(1:iscn-1));
            tmx( tmx==ssg+nTc | tmx==esg+nTc)=[]; 
        end
        nmx=length(tmx);
        
        cn=c-sum(c)/nT; cn=cn/sqrt(cn*cn');
        if cn*cn1>e || cn*cn2>e || cn*cn3>e, break; end
        cn3=cn2; cn2=cn1; cn1=cn'; itr=itr+1; 
    end
%     clear cn cn1 cn2 cn3
    
    if nmx>1
        T=zeros(nX,PLe,'single');
        tS=tmx-PLh(1); tE=tmx+PL-1+PLh(2); 
        for i=1:nmx, ts=tS(i); te=tE(i);
            zs=[]; if ts<=0, zs=zeros(nX,abs(ts)+1,'single'); ts=1; end
            ze=[]; if te>nT, ze=zeros(nX,te-nT,'single'); te=nT; end
            T=T+[zs D(:,ts:te) ze];
        end; TMPL{itp}=T/nmx;
%         clear T
        TMXTMPL{itp}=single(tmx); CTMPL(itp,:)=single(c); 
        SCMX(itp)=single(sum(c(tmx))); ITR(itp)=itr;
    end
    if ~rem(itp,stp), fprintf(['%s ' s '-ITP%d\n'], datestr(now,'HH:MM:SS'), itp); end
end
% if rem(itp,stp), fprintf([s '-ITP%d\n'],itp); end
clear SGf SGfn c
 
%% Selecting QPP
[~,ITP1]=max(SCMX); 
QPP=TMPL{ITP1}; TMX=TMXTMPL{ITP1}; C=CTMPL(ITP1,:); ITER=ITR(ITP1);
MET=[median(C(TMX)) median(diff(TMX))*tres length(TMX)];



