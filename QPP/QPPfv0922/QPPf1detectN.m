function [QPPn, TMXn, METn]=QPPf1detectN(D, C, ntlist, PL, cth,tres)
PLh=round(PL/2)+[0 -rem(PL,2)]; % pad length for temporally extending a QPP
% PLc=round(PL/2)+(1:PL); % range of interest in an extended QPP, matches PLh
% PLe=PL+sum(PLh); % length of an extended QPP (derivable but saved/read)
% tsh=floor(PL/4); % max timeshift when comparing QPPs or any 2 templates
ssg=1;  % starting segment/scan ~, SEE QPPf4regscn
    
c=double(C); nT=size(D,2); [nX,nT]=size(D); nXL=nX*PL; PLe=PL+sum(PLh);
nscn=length(ntlist); ith=2;
[~,tmx]=findpeaks(-c,'MinPeakHeight',cth(ith),'MinPeakDistance',PL);
for iscn=1:nscn
    nt=ntlist(iscn); esg=nt-PL+1; nTc=sum(ntlist(1:iscn-1));
    tmx( tmx==ssg+nTc | tmx==esg+nTc)=[]; 
end; nmx=length(tmx);

T=zeros(nX,PLe,'single');
tS=tmx-PLh(1); tE=tmx+PL-1+PLh(2); 
for i=1:nmx, ts=tS(i); te=tE(i);
    zs=[]; if ts<=0, zs=zeros(nX,abs(ts)+1,'single'); ts=1; end
    ze=[]; if te>nT, ze=zeros(nX,te-nT,'single'); te=nT; end
    T=T+[zs D(:,ts:te) ze];
end

QPPn=T/nmx; TMXn=tmx; METn=[median(C(TMXn)) median(diff(TMXn))*tres length(TMXn)];

    