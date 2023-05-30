function [RPs,RPs2,RPs3]=sig_str(RPp)
    TSIG1=0.05;
    TSIG2=0.01;
    TSIG3=0.001;
    RPs=' ';
    if RPp<TSIG3
        RPs='***';
        RPs2=RPs;
    elseif RPp<TSIG2
        RPs='**';
        RPs2=RPs;
    elseif RPp<TSIG1
        RPs='*';
        RPs2=RPs;
    else
        RPs2='n.s.';
        RPs=' ';
    end
    if RPp>1
        RPs3=sprintf('%s (>1)',RPs2);
    elseif RPp<0.0001
        RPs3=sprintf('%s (<0.0001)',RPs2);
    else
        RPs3=sprintf('%s (%3.1g)',RPs2,RPp);
    end
    
end