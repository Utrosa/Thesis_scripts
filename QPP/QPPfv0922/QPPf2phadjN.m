function [QPPa,TMXa,Ca,METa,SDa,SD,flga,cT1Tj,nsim]=...
    QPPf2phadjN(QPP,TMPL,SCMX,TMXTMPL,CTMPL,tsh,PLc,cthph,sdph,s,tres)
%%%% added by Nan Xu for detecting the phase adjusted reverse QPP template.
nX=size(QPP,1); PL=length(PLc); nXL=nX*PL; 
nITP=length(TMPL); nsd=length(sdph); nT=size(CTMPL,2); PLe=size(QPP,2);

%% Comparing QPP with other templates, finding similar ones & sorting
T1=Tshcfn(QPP,tsh,PLc); % shifting QPP, flattening & normalizing
Tj=zeros(nXL,nITP,'single'); % flattening & normalizing other templates
parfor i=1:nITP, T=TMPL{i}; if any(T), T=T(:,PLc); T=T(:); T=T-sum(T)/nXL; 
        Tj(:,i)=T'/sqrt(T'*T); end; end
cT1Tj=max(abs(T1*Tj))'; % max corr btwn shifted QPP & other templates
clear T1 Tj

ITPsim=find(cT1Tj>=cthph); % similar templates to QPP
[~,isscmx]=sort(SCMX(ITPsim),'descend'); 
I=ITPsim(isscmx); % similar templates sorted based on sum C@maxima
nsim=length(ITPsim); 

%% Finding seed timecourses (y) & applying strict(s) & relaxed(r) criteria
z=zeros(1,nsd,'single'); ITPs=z; ITPr=z; flga=z;
z=cell(1,nsd); z(:)={single([])}; QPPa=z; TMXa=z;
Ca=zeros(nsd,nT,'single'); METa=zeros(nsd,3,'single'); 
z=zeros(nsd,PLe,'single'); SDa=z; SD=z; clear z

for isd=1:nsd
    y=zeros(nsim,PL,'single'); t=zeros(nsim,2,'single');  
    for i=1:nsim
        y(i,:)=mean(TMPL{I(i)}(sdph{isd},PLc),1);
        [~,t(i,1)]=max(y(i,:)); [~,t(i,2)]=min(y(i,:));
    end
    
    if s
        for i=1:nsim
            if abs(y(i,end))<0.1 && mean(y(i,end-3:end))>0.1 && ...
                    t(i,1)>=PL/2 && t(i,1)>t(i,2)
                ITPs(isd)=I(i); break
            end
        end
        if ~ITPs(isd)
            for i=1:nsim
                if abs(y(i,end))<0.2 && mean(y(i,end-3:end))>0.2 && ...
                        t(i,1)>=PL/2 && t(i,1)>t(i,2)
                    ITPs(isd)=I(i); break
                end
            end
        end
    end

    ITPr(isd)=ITPs(isd);
    if ~ITPr(isd)
        for i=1:nsim
            if t(i,1)>=PL/2 && t(i,1)>t(i,2), ITPr(isd)=I(i); break; end
        end
    end

    if ITPr(isd)
        QPPa{isd}=TMPL{ITPr(isd)}; 
        TMXa{isd}=TMXTMPL{ITPr(isd)}; 
        Ca(isd,:)=CTMPL(ITPr(isd),:);
        METa(isd,:)=[median(Ca(isd,TMXa{isd})) ...
            median(diff(TMXa{isd}))*tres length(TMXa{isd})]; 
        SDa(isd,:)=mean(QPPa{isd}(sdph{isd},:),1); 
        flga(isd)=any(ITPs(isd));
    end; SD(isd,:)=mean(QPP(sdph{isd},:),1);
end

