function A=notnan_percentile(B,P,percent)
assert(percent>=0)
assert(percent<=1)

assert(P==1); %not important but identical in format to mean(B,1)
A=nan(1,size(B,2));
for I=1:size(B,2)
    vec=B(:,I);
    vec=vec(~isnan(vec));
    vec=vec(~isinf(vec));
    if isempty(vec)
        vec=nan;
    end
    
    pos=min(floor(percent*length(vec))+1,length(vec));
    svec=sort(vec);
    
    A(1,I)=svec(pos);
end
