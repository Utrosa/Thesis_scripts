function [matchesD,peaks_num,BOOTS]=boot_and_match_best(mydata,BOOTS_REPO,LOCSD,PKSD,xx,KernelWidth,TH_PROX)

NBOOT=size(BOOTS_REPO,1);
matchesD=nan(NBOOT,length(LOCSD));
peaks_num=nan(NBOOT,1);

BOOTS=cell(NBOOT,1);

IS_PLOT=false; % never plot!
for B=1:NBOOT
    if mod(B,10)==1
        fprintf('.')
    end
    boot=BOOTS_REPO{B}.boot;
    datB=mydata(boot,:,:);

    HB=compute_marginal(datB,xx,KernelWidth);
    [PKSB,LOCSB]=findpeaks(HB,xx, 'Annotate','peaks');

    if IS_PLOT % never happens (for debugging)
        plot(xx,HB,'-r', 'LineWidth',0.5);hold on;
        plot(LOCSB,PKSB,'x');hold on;
        plot(LOCSB,PKSB,'ro','MarkerFaceColor','r');hold on;
    end

    matchB=nan(size(LOCSB,1),2);
    for ll=1:length(LOCSB)
        loc=LOCSB(ll);
        loc_pk=PKSB(ll);
        candidats=LOCSD(abs(LOCSD-loc)<TH_PROX);
        [~,idx]=min(abs(LOCSD-loc));
        candidat=LOCSD(idx);
        candidat_pk=PKSD(idx);
        assert(length(candidat)==1);
        if abs(candidat-loc)<TH_PROX
            matchB(ll,:)=[candidat,candidat_pk];
            if isnan(matchesD(B,idx))
                matchesD(B,idx)=loc;
            else
                if abs(candidat-loc)<abs(matchesD(B,idx)-candidat)
                    matchesD(B,idx)=loc;
                end
            end


            if IS_PLOT % never happens (for debugging)
                plot([loc,candidat],[loc_pk,candidat_pk],'g-','LineWidth',3);hold on;
                plot(loc,loc_pk,'gd','MarkerFaceColor','g');hold on;
            end

        end
    end
    peaks_num(B)=length(LOCSB);
    BOOTS{B}.match_peak=matchB(:,1);
    BOOTS{B}.HB=HB;
    BOOTS{B}.PKSB=PKSB;
    BOOTS{B}.LOCSB=LOCSB;
    BOOTS{B}.peak_shared_origins=LOCSD;
    BOOTS{B}.peak_num=length(LOCSB);

end
fprintf('\n')