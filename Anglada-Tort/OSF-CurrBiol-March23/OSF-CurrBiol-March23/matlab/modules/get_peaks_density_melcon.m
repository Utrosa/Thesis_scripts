function [pks_density,PKSD0,LOCSD0]=get_peaks_density_melcon(xvals,yvals,BOOTS_REPO,xx,KernelWidth)
fprintf('get peaks density...')
pks_density=zeros(size(xx));
NBOOT=size(BOOTS_REPO,1);
for B=1:NBOOT
    if mod(B,10)==1
        fprintf('.')
    end
    boot=BOOTS_REPO{B}.boot;

    [HB,~]=smooth_dense_data(xvals(boot),yvals(boot),xx,KernelWidth);


    %HB=compute_marginal(datB,xx,KernelWidth);
    [~,LOCSB]=findpeaks(HB,xx, 'Annotate','peaks');

    H=zeros(size(xx));
    for ll=1:length(LOCSB)
        h=normpdf(xx,LOCSB(ll),KernelWidth);
        h=h/sum(h);
        H=H+h;
    end

    pks_density=pks_density+H;

end
pks_density=pks_density/sum(pks_density);
[PKSD0,LOCSD0]=findpeaks(pks_density,xx, 'Annotate','peaks');
fprintf('\n')