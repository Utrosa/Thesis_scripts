function [pyx,pyv,yvar,yvar2,p1,yp,yvv,yvv2]=compute_biases(xvals,yvals,xx,KW,VAR_TH_RANGE,MODEL_DEG)

pos=~isnan(xvals+yvals);
xvals=xvals(pos);
yvals=yvals(pos);


pyx=nan(size(xx));
pyv=nan(size(xx));

for xxx=1:length(xx)
    x=xx(xxx);
    w=exp(- ((x-xvals).^2 )/ (KW^2) );
    w=w/sum(w);
    interp_val=sum(w.*yvals);
    interp_var=sqrt(sum(w.*(yvals-interp_val).^2));
    pyx(xxx)=interp_val;
    pyv(xxx)=interp_var;

end
yvar=std(yvals(:));
yvar2=std(yvals(abs(xvals)<VAR_TH_RANGE));


p1 = polyfit(xvals,yvals,MODEL_DEG);
yp = polyval(p1,xx);


yvv=ones(size(xx))*yvar;
yvv2=ones(size(xx))*yvar2;

