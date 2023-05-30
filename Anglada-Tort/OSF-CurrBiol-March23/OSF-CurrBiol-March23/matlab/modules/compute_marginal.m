function H=compute_marginal(dat,xx,KernelWidth)
H=zeros(size(xx));
marg=dat(:,:,:);
marg=marg(:);
marg=marg(~isnan(marg));
for ll=1:length(marg)
    h=normpdf(xx,marg(ll),KernelWidth);
    h=h/sum(h);
    H=H+h;
end

H=H/sum(H);

