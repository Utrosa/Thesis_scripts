function paramQPPf1=param_QPPf1(nP, ntlist, PL, cth13, cth45, ssg)
    cth=cell(nP,1); cth(1:min(nP,3))={cth13}; cth(4:end)={cth45}; % correlation threshold in the main algorithm
    ncth1=3*ones(nP,1); if nP>=4, ncth1(4:end)=3; end % #iters with lower cth
    nitr=15; % max #iterations in the main algorithm

    ITP=cell(nP,1); ITP1=cell(nP,1); PLh=cell(nP,1);
    t=zeros(1,nP); t(1:min(nP,3))=8; t(4:end)=5; % for fast detection of QPP
    for ip=1:nP
        [ITP{ip}, ITP1{ip}]=ITPgenerate(ntlist, PL(ip), ssg(ip), t(ip));
        [PLh{ip}, ~,~]=PLextension(PL(ip));
    end
    
    paramQPPf1.PL=PL; 
    paramQPPf1.cth=cth; paramQPPf1.ncth1=ncth1; paramQPPf1.nitr=nitr;
    paramQPPf1.ssg=ssg; paramQPPf1.ITProbust=ITP; paramQPPf1.ITPfast=ITP1; paramQPPf1.PLh=PLh;
