function d=JSD2(P,Q)
    assert(numel(P)==numel(Q));
    assert(sum(size(P)==size(Q))==length(size(P)));
    P=P(:);
    Q=Q(:);
    P=P/sum(P);
    Q=Q/sum(Q);
    M=(P+Q)/2;
    d= (DKL2(P,M)+DKL2(Q,M))/2;
    

    