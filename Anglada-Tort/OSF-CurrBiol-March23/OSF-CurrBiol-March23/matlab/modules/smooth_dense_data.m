function [pyx,pyv]=smooth_dense_data(xvals,yvals,xx,KW)
pyx=nan(size(xx));
pyv=nan(size(xx));
for ll=1:length(xx)

%     if mod(ll,round(length(xx)/10))==1
%         fprintf('.')
%     end

    x=xx(ll);
    w=exp(- ((x-xvals).^2 )/ (2*(KW^2)) );
    w=w/sum(w);
    interp_val=sum(w.*yvals);
    interp_var=sqrt(sum(w.*(yvals-interp_val).^2));
    pyx(ll)=interp_val;
    pyv(ll)=interp_var;
end