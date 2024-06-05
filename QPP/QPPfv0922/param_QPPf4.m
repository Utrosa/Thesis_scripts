function paramQPPf4=param_QPPf4(nP, PL, ROI2Net, fz)    
    PLc=cell(nP,1); 
    for ip=1:nP
        [~, PLc{ip}, ~]=PLextension(PL(ip));
    end
    
    nnet=length(unique(ROI2Net)); iROI2NetC=cell(nnet,1); iNetL=zeros(nnet+1,1);
    [a,iROI2Net]=sort(ROI2Net); 
    for inet=1:nnet
        iROI2NetC{inet}=iROI2Net(a==inet);
        iNetL(inet+1)=iNetL(inet)+length(iROI2NetC{inet});
    end
    paramQPPf4.PLc=PLc;
    paramQPPf4.iNetL=iNetL;
    paramQPPf4.iROI2NetC=iROI2NetC;
    paramQPPf4.fz=fz;
    
%     paramQPPf4.nnet=nnet;
    paramQPPf4.iROI2Net=iROI2Net;
    
