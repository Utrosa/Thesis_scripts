function paramQPPf2=param_QPPf2(nP, PL, cthph, sdph, s)      
    PLc=cell(nP,1); tsh=zeros(nP,1);
    for ip=1:nP
        tsh(ip)=floor(PL(ip)/4); % max timeshift when comparing 2 QPP templates  
        [~, PLc{ip}, ~]=PLextension(PL(ip));
    end
    nsd=zeros(nP,1); for ip=1:nP, nsd(ip)=length(sdph{ip}); end; nSD=max(nsd);
    
    paramQPPf2.tsh=tsh; 
    paramQPPf2.PLc=PLc;
    paramQPPf2.cthph=cthph;
    paramQPPf2.sdph=sdph;
    paramQPPf2.s=s;
    paramQPPf2.nSD=nSD;
