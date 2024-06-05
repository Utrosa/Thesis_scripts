function To=Tshcfn(T,tsh,PLc)
nX=size(T,1); nXL=nX*length(PLc);
To=zeros(2*tsh+1,nXL,'single'); 
To(tsh+1,:)=Tcfn(T,PLc);
for i=1:tsh
    To(i,:)=Tcfn([T(:,tsh+1-i+1:end) zeros(nX,tsh+1-i,'single')],PLc);
    To(tsh+1+i,:)=Tcfn([zeros(nX,i,'single') T(:,1:end-i)],PLc);
end