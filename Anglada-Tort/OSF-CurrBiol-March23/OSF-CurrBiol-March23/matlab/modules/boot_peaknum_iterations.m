function [peaks_num_i,BOOTSi]=boot_peaknum_iterations(mydata,BOOTS_REPO,xx,KernelWidth,LIMIT_RANGE)
NBOOT=size(BOOTS_REPO,1);
NK=size(mydata,2);
BOOTSi=cell(NBOOT,1);
peaks_num_i=nan(NBOOT,NK);
for B=1:NBOOT
    if mod(B,10)==1
        fprintf('.')
    end
    boot=BOOTS_REPO{B}.boot;

    for KK=1:NK
        datBi=mydata(boot,KK,:); %from all data!

        HBi=compute_marginal(datBi,xx,KernelWidth);

        [PKSBi,LOCSBi]=findpeaks(HBi,xx, 'Annotate','peaks');
        if ~isnan(LIMIT_RANGE)
            pos=(LOCSBi<=max(LIMIT_RANGE))&(LOCSBi>=min(LIMIT_RANGE));
            PKSBi=PKSBi(pos);
            LOCSBi=LOCSBi(pos);
        end

        peaks_num_i(B,KK)=length(LOCSBi);
    end

    BOOTSi{B}.HBi=HBi;
    BOOTSi{B}.PKSBi=PKSBi;
    BOOTSi{B}.LOCSBi=LOCSBi;
    BOOTSi{B}.peak_num=length(LOCSBi);

end
fprintf('\n')