function d=DKL2(P,Q)
    assert(abs(sum(P(:))-1)<1e-12);
    assert(abs(sum(Q(:))-1)<1e-12);
    
    d=sum(P(:).*log2((P(:)+eps)./(Q(:)+eps)));
    
    