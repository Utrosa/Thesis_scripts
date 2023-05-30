function plot_simulations(HASi,HAS_sim,xxx,XRANGE,mname)

NK1=length(HASi);

vmax=-9999;
for K=1:NK1
    vmax=max(vmax,HASi{K}(2:(end-1)));
    vmax=max(vmax,HAS_sim{K}(2:(end-1)));
end
vmax=max(max(vmax));


myKS=1:2:NK1;


NSUB1=2;
NSUB2=6;
legs={'singing data',mname};
for KK=1:length(myKS)

    K=myKS(KK);
    subplot(NSUB1,NSUB2,KK)
    plot(xxx,HASi{K},'k','LineWidth',4);
    hold all;title(sprintf('%s - iter=%d',legs{1},K-1));
    
    set(gca,'Xtick',min(XRANGE):max(XRANGE))
    xlim([min(XRANGE)+0.3,max(XRANGE)-0.3]);
    ylim([0 vmax]);
    set(gca,'FontSize',14)
    mscore=JSD2(HASi{K}(2:(end-1))',HAS_sim{K}(2:(end-1)));


    subplot(NSUB1,NSUB2,KK+1*NSUB2)
    plot(xxx,HAS_sim{K},'r','LineWidth',2);hold all;
    hold all;title(sprintf('%s - %3.3f ',legs{2},mscore));
    set(gca,'Xtick',min(XRANGE):max(XRANGE))
    xlim([min(XRANGE)+0.3,max(XRANGE)-0.3]);
    ylim([0 vmax]);
    set(gca,'FontSize',14)


end