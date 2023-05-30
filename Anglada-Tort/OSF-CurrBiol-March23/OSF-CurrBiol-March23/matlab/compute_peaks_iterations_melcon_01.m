addpath('support-functions/'); % dependency on internal libraries
addpath('modules/'); % dependency on modules

close all;clc;clear all;
tic()

output_dir='output/';
mkdir(output_dir);
IS_GET_PEAK_DENSITY_FOR_PEAKS=true;

EXPRS={'melcon-2tones'};

KernelWidth=0.25;
TH_PROX=0.5;
PEAK_SIG_TH=[0.95,0.90,0.50]; %thershold for signficance
%STD_FACT=1;
STD_FACT=1.96; % param for display

xx0=-20:0.1:20;
INT_RANGE=-16:16;MM=5; % params for integer score
NBOOT=1000;
IS_PLOT=true;
IS_SAVE=true;
XRANGE=[-15,15];
RES=cell(size(EXPRS));
LIMIT_RANGE=[-20,20]; % limit the counting of peaks to this range

for EE=1:length(EXPRS)
    EXPR=EXPRS{EE};

    xx=xx0;
    NI_MAX=99999;

    DAT=[];

    [data,raw_data,UN,NK,NI,unets]=analyze_curbio_06(EXPR);

    for III=0:NI

        RANGE=1:NI;
        RANGE_msg='all';
        fprintf('Now in %s %s\n---------------------------\n',EXPR,RANGE_msg);
        %%

        xvals=data(:,1,1);
        yvals=data(:,2,1);
        pos=~isnan(xvals+yvals);
        xvals=xvals(pos);
        yvals=yvals(pos);
        pos2=abs(xvals)<max(XRANGE);
        xvals=xvals(pos2);
        yvals=yvals(pos2);

        KW=KernelWidth;
        [HA,HAv]=smooth_dense_data(xvals,yvals,xx,KW);
        [PKSA,LOCSA]=findpeaks(HA,xx, 'Annotate','peaks');

        if IS_PLOT
            figure(100+EE);clf


            plot(xx,HA,'k-','Linewidth',2);
            hold on;
            plot(LOCSA,PKSA,'ob','MarkerFaceColor','r');
            for ll=INT_RANGE
                plot([ll ll],[min(HA),max(HA)],'c-');hold on;

            end
        end

        BOOTS_REPO=cell(NBOOT,1);
        fprintf('randomizing once bootrapps %s %s\n',EXPR,RANGE_msg)
        %randomized once  bootrapping
        for B=1:NBOOT
            boot=randi(size(xvals,1),size(xvals,1),1);
            BOOTS_REPO{B}.boot=boot;
        end

        pks_density=[];

        if IS_GET_PEAK_DENSITY_FOR_PEAKS
            [pks_density,PKSD0,LOCSD0]=get_peaks_density_melcon(xvals,yvals,BOOTS_REPO,xx,KernelWidth);
            PKSD_interp=interp1(xx,HA,LOCSD0);

            LOCSD=LOCSD0;
            PKSD=PKSD_interp;

        else
            LOCSD=LOCSA;
            PKSD=PKSA;

        end
        IS_PLOT2=false;
        fprintf('Bootstrapping peaks [%s %s] ...\n',EXPR,RANGE_msg)
        [matchesD,peaks_num,BOOTS]=boot_and_match_best_melcon(xvals,yvals,BOOTS_REPO,LOCSD,PKSD,xx,KernelWidth,TH_PROX,IS_PLOT2);


        fprintf('computing peaks stats [%s %s] ...\n',EXPR,RANGE_msg)
        PEAK_SIG=sum(~isnan(matchesD))/NBOOT;
        PEAKS_MEAN=notnan_mean(matchesD,1);
        PEAKS_STD=STD_FACT*notnan_std(matchesD,1);
        PEAKS_STD(PEAKS_STD<1e-10)=nan;
        PEAKS_INT=round(PEAKS_MEAN);
        PEAKS_Z=abs(PEAKS_MEAN-PEAKS_INT)./(PEAKS_STD+eps);
        PEAKS_D=abs(PEAKS_MEAN-PEAKS_INT);


        figure(100+EE);
        plot_summary_graph(EXPR,xx,HA,INT_RANGE,PEAKS_MEAN,PEAKS_STD,PEAK_SIG,PEAKS_Z,PEAKS_D,PEAK_SIG_TH,LOCSD,PKSD,LOCSA,PKSA);
        title(sprintf ("%s %s",EXPR,RANGE_msg));
        ylim([0.9*min(HA),1.1*max(HA)]);

        fprintf('saving final data [%s %s] ...\n',EXPR,RANGE_msg)

        DAT.data=data;
        DAT.raw_data=raw_data;
        DAT.UN=UN;
        DAT.NK=NK;
        DAT.NI=NI;
        DAT.unets=unets;
        DAT.BOOTS=BOOTS;

        DAT.BOOTS_REPO=BOOTS_REPO;
        DAT.NBOOT=NBOOT;
        DAT.pks_density=pks_density;
        DAT.PKSD=PKSD; % saved possible locations of peaks
        DAT.LOCSD=LOCSD;% saved possible locations of peaks
        DAT.PKSA=PKSA; % all data
        DAT.LOCSA=LOCSA;% all data
        DAT.peaks.PEAK_SIG=PEAK_SIG;
        DAT.peaks.PEAKS_MEAN=PEAKS_MEAN;
        DAT.peaks.PEAKS_STD=PEAKS_INT;
        DAT.peaks.PEAKS_Z=PEAKS_Z;
        DAT.peaks.PEAKS_D=PEAKS_D;
        DAT.peaks.origin_location= LOCSD;
        DAT.peaks.origin_estimated_density= PKSD;
        DAT.peaks.matched_locations=matchesD;
        DAT.xx=xx;


        RES{EE}.DAT=DAT;

        if IS_SAVE

            obfname=sprintf('%s/boottsrapping_peaks_matrix-%s-%s.NBOOT-%d.txt',output_dir,EXPR,RANGE_msg,NBOOT);
            fprintf('saving bootrapping data to filename: %s\n',obfname)
            writematrix(matchesD,obfname)

            ojfname=sprintf('%s/boottsrapping_peaks_info-%s-%s.NBOOT-%d.json',output_dir,EXPR,RANGE_msg,NBOOT);
            fprintf('saving peaks info to a json: %s\n',ojfname)
            js=jsonencode(DAT.peaks);
            FID=fopen(ojfname,"w");
            fprintf(FID,'%s',js);
            fclose(FID);


        end

        fprintf('finished for %s %s!!\n\n',EXPR,EXPR,RANGE_msg);
    end

    if IS_PLOT
        figure(100+EE);
        set(gcf,'Units','normalized');
        set(gcf,'Position',[0 0  1.0000    0.9143]);
        drawnow;

        if IS_SAVE
            fprintf('saving exported data... %s \n',EXPR)
            figure(100+EE);
            ofname=sprintf('%s/summary_plot-%s.NBOOT-%d.png',output_dir,EXPR,NBOOT);
            fprintf('saving png image to filename: %s\n',ofname)
            print(ofname, '-dpng'); % save png file

        end

    end
end

toc()
