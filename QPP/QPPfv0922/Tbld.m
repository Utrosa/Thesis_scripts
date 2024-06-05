function T=Tbld(D,tmx,PL,PLh,dv)
[nX,nT]=size(D); PLe=PL+sum(PLh); nmx=length(tmx); 
tS=tmx-PLh(1); tE=tmx+PL-1+PLh(2);
T=zeros(nX,PLe,'single'); 
for i=1:nmx
    ts=tS(i); te=tE(i);
    zs=[]; if ts<=0, zs=zeros(nX,abs(ts)+1,'single'); ts=1; end
    ze=[]; if te>nT, ze=zeros(nX,te-nT,'single'); te=nT; end
    T=T+[zs D(:,ts:te) ze];
end
if dv, T=T/nmx; end