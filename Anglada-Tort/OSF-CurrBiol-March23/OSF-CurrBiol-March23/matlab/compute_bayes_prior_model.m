function [pRR,pPkS]=compute_bayes_prior_model(pR,sigma_sen,p0,xxx,NK1)
M=length(p0);
assert(length(pR)==M);

pSR=nan(M,M);
    for I=1:M
        z=xxx(I);
        Q=normpdf(xxx,z,sigma_sen);Q=Q./sum(Q(:));
        pSR(:,I)=Q;
    end


%%% perform Bayesian inference p(R^|S) for this chain R->S->R^
pRS=nan(M,M);
for I=1:M
    pRS(:,I)=pSR(I,:).*(pR');
    pRS(:,I)=pRS(:,I)/sum(pRS(:,I));
end

pRR=pRS*pSR; % One step of the iterative process stimuli--> noise --> inference (P(Rn+1|Rn) for the iterative experiment)

pPkS=cell(NK1,1);

pPk=p0;K=1;
pPkS{K}=pPk;
for K=2:NK1
    pPk=pRR*pPk;
    pPkS{K}=pPk;
end
