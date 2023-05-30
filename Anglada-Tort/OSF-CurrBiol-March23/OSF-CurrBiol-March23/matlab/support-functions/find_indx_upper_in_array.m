function idx=find_indx_upper_in_array(vals,xx)
% find the smallest value in xx that is larger than vals

N=length(xx);
M=length(vals);
assert(prod(size(vals))==length(vals)); % a vector;
assert(prod(size(xx))==length(xx)); % a vector;
xx=reshape(xx,N,1);
vals=reshape(vals,1,M);

MAT=-(repmat(vals,N,1)-repmat(xx,1,M)); 
MAT(MAT<0)=nan;% only places where xx is larger than vals

[~,idx]=min(MAT);
