function To=Tcfn(T,PLc)
T=T(:,PLc); T=T(:); T=T-sum(T)/length(T); To=T'/sqrt(T'*T);