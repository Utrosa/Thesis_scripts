function idx=find_indx_in_array(vals,xx)


N=length(xx);
M=length(vals);
assert(prod(size(vals))==length(vals)); % a vector;
assert(prod(size(xx))==length(xx)); % a vector;
xx=reshape(xx,N,1);
vals=reshape(vals,1,M);


[~,idx]=min(abs(repmat(vals,N,1)-repmat(xx,1,M)));
