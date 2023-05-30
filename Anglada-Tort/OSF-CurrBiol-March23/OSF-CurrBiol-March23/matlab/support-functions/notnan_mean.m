function A=notnan_mean(B,P)
assert(P==1); %not important but identical in format to mean(B,1)
A=nan(1,size(B,2));
for I=1:size(B,2)
    vec=B(:,I);
    vec=vec(~isnan(vec));
    vec=vec(~isinf(vec));
    if isempty(vec)
        vec=nan;
    end
    A(1,I)=mean(vec);
end
