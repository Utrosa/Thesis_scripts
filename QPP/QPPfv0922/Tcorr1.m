function [c,tmx]=Tcorr1(D,ntlist,T,Tfn,nscn,ssg,cth)
[nX,nT]=size(D); nXL=max(length(T(:)),length(Tfn)); PL=nXL/nX; 

if any(T), T=T(:); T=T-sum(T)/nXL; Tfn=T'/sqrt(T'*T); end

c=zeros(nT,1);
for i=1:nscn
    nt=ntlist(i); esg=nt-PL+1; nT0=sum(ntlist(1:i-1));
    for isg=ssg:esg
        t=isg+nT0; S=D(:,t:t+PL-1);        
        S=S(:); S=S-sum(S)/nXL; S=S/sqrt(S'*S); c(t)=Tfn*S;
    end
end

tmx=[];
if cth
    [~,tmx]=findpeaks(c,'MinPeakHeight',cth,'MinPeakDistance',PL);
    for i=1:nscn, tmx( tmx==ssg+(i-1)*nt | tmx==esg+(i-1)*nt)=[]; end
end

c=single(c); tmx=single(tmx);
