function [QPP,C,mC,QPPa,Ca,mCa]=QPPf3xresidRbst(D, ntlist, TMX,TMXa,PLh,PLc,nscn,ssg)
PL=length(PLc); nsd=length(TMXa); nT=size(D,2); 

QPP=Tbld(D,TMX,PL,PLh,1); 
C=Tcorr1(D,ntlist,QPP(:,PLc),[],nscn,ssg,0);
mC=median(C(TMX));

QPPa=cell(1,nsd); QPPa(:)={single([])}; 
Ca=zeros(nsd,nT,'single'); mCa=zeros(nsd,1,'single');
for i=1:nsd
    tmx=TMXa{i};
    if any(tmx)
        QPPa{i}=Tbld(D,tmx,PL,PLh,1); 
        Ca(i,:)=Tcorr1(D,ntlist,QPPa{i}(:,PLc),[],nscn,ssg,0);
        mCa(i,1)=median(Ca(i,tmx));
    end
end