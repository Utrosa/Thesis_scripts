function [c,tmx]=Tcorr(D,T,Tfn,nscn,ssg,cth)
[nX,nT]=size(D); nXL=max(length(T(:)),length(Tfn)); PL=nXL/nX; 
nt=nT/nscn; esg=nt-PL+1;

if any(T), T=T(:); T=T-sum(T)/nXL; Tfn=T'/sqrt(T'*T); end

c=zeros(nT,1);
for i=1:nscn
for isg=ssg:esg
    t=isg+(i-1)*nt; S=D(:,t:t+PL-1);
    S=S(:); S=S-sum(S)/nXL; S=S/sqrt(S'*S); c(t)=Tfn*S;
end
end

tmx=[];
if cth
    [~,tmx]=findpeaks(c,'MinPeakHeight',cth,'MinPeakDistance',PL);
    for i=1:nscn, tmx( tmx==ssg+(i-1)*nt | tmx==esg+(i-1)*nt)=[]; end
end

c=single(c); tmx=single(tmx);
