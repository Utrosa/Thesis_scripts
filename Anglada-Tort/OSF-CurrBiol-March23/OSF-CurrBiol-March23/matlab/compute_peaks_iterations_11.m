addpath('support-functions/'); % dependency on internal libraries
addpath('modules/'); % dependency on modules

close all;clc;clear all; % "To begin at the beginning"/ Under Milk Wood – Dylan Thomas

tic() % count time of the entire process

output_dir='output/'; % directory for output
mkdir(output_dir);


IS_GET_PEAK_DENSITY_FOR_PEAKS=true; % this means we use pre-computations of peak density 


EXPRS={'btw-3tones','btw-2tones', 'wth-3tones','btw-5tones','btw-free','slider-2tones','mem-3tones-10','mem-3tones-5','mem-3tones-match','btw-3tones-india','wth-3tones-india'}; % list of all experiments to process



%%% Paramters of the experiment
KernelWidth=0.25; 
TH_PROX=0.5; % maximal distance in semitone between a given peak and it's associated category
NBOOT=1000; % number of bootrapping datasets
xx0=-20:0.01:20; % the range of intervals to consider the step (0.01) is the resolution for most computations.
INT_RANGE=-16:16; % integer range to work with
IS_SAVE=true; % if false do not save output
LIMIT_RANGE=[-20,20]; % limit the counting of peaks to this range 


% display parameters:
IS_PLOT=true; % plot on screen?
PEAK_SIG_TH=[0.95,0.90,0.50]; %thershold for signficance in terms of plotting the peaks
STD_FACT=1.96; % param for display - matching STD in terms of CI.

RES=cell(size(EXPRS)); % container for all experimental results

parfor EE=1:length(EXPRS) % use parallel processing to get all data
    EXPR=EXPRS{EE};

    %%% only for btw-free expriment we compute just the first 3 marginals
    %%% (there are just too many for this experiment) in all other casees
    %%% we compute all marginals.
    if strcmp(EXPR,'btw-free')
        xx=-32:0.01:32;
        NI_MAX=3; 
    else
        xx=xx0;
        NI_MAX=99999;
    end    


    DAT=[];% container for the data for this iteration

    [data,raw_data,UN,NK,NI,unets]=analyze_curbio_06(EXPR); % read the data from the file.
    if IS_PLOT
        figure(100+EE);clf
    end

    % set subplots for dispaly - choose one that allow you to show nicely
    % aall iterations.
    NSUB1= 2;
    NSUB2=floor((NI+2.01)/NSUB1);
    if NSUB2*NSUB1<(NI+2)
        NSUB2=NSUB2+1;
    end
    if strcmp(EXPR,'btw-free')
        NSUB1= 2;
        NSUB2= 3;
    end

    % go over all intervals:
    for III=0:NI
        if III>NI_MAX % in case of the free experiment skip all iterations except the first 3.
            fprintf("for %s experiment I skip interval marginal larger than 3 (skipping %d)",EXPR,III)
            continue
        end
        
        % we want to comptue marginals for all intervals, but ALSO compute
        % the marginals for all intervals together (this is conidition III==0)
        if III==0 % in this case use ALL intervals
            RANGE=1:NI;
            RANGE_msg='all-intervals'; % choose all intervals
        else
            RANGE=III; % choose only one interval
            RANGE_msg=sprintf('interval-%d',III);
        end

        fprintf('Now in %s %s\n---------------------------\n',EXPR,RANGE_msg); % progress report

        dat=data(:,(NK-2):NK,RANGE); % this line take data from the last 3 iterations (NK-2):NK) and "RANGE" which is the intervals that we want (RANGE can be all intervals).
    
        HA=compute_marginal(dat,xx,KernelWidth); % we compute the marginal distribution
        [PKSA,LOCSA]=findpeaks(HA,xx, 'Annotate','peaks'); % we find peaks for the overall marginal

        BOOTS_REPO=cell(NBOOT,1); % we randomized once bootrapping data - we will reuse this bootrapping datsest through out.
        
        
        %randomized once  bootrapping
        fprintf('randomizing once bootrapps %s %s\n',EXPR,RANGE_msg)
        for B=1:NBOOT
            boot=randi(size(data,1),size(data,1),1);
            BOOTS_REPO{B}.boot=boot;
        end

        
        % gets the peak density from the bootrapping datasets. based on it
        % derive the candidate cateogry locations
        pks_density=[];
        if IS_GET_PEAK_DENSITY_FOR_PEAKS
            [pks_density,PKSD0,LOCSD0]=get_peaks_density(dat,BOOTS_REPO,xx,KernelWidth);
            PKSD_interp=interp1(xx,HA,LOCSD0);

            LOCSD=LOCSD0;
            PKSD=PKSD_interp;

        else
            LOCSD=LOCSA;
            PKSD=PKSA;

        end
        
        % find for each peak the category location associated with it
        fprintf('Bootstrapping peaks [%s %s] ...\n',EXPR,RANGE_msg)
        %OLD: [matchesD,peaks_num,BOOTS]=boot_and_match(dat,BOOTS_REPO,LOCSD,PKSD,xx,KernelWidth,TH_PROX);
        [matchesD,peaks_num,BOOTS]=boot_and_match_best(dat,BOOTS_REPO,LOCSD,PKSD,xx,KernelWidth,TH_PROX);

        % bootrapping iteration data to get the number of peaks per
        % iterations
        fprintf('Bootstrapping peaks iterations  [%s %s] ...\n',EXPR,RANGE_msg)
        [peaks_num_i,BOOTSi]=boot_peaknum_iterations(data,BOOTS_REPO,xx,KernelWidth,LIMIT_RANGE);

        % compute category peaks statistics:    
        fprintf('computing peaks stats [%s %s] ...\n',EXPR,RANGE_msg)
        PEAK_SIG=sum(~isnan(matchesD))/NBOOT;
        PEAKS_MEAN=notnan_mean(matchesD,1);
        PEAKS_STD=STD_FACT*notnan_std(matchesD,1);
        PEAKS_STD(PEAKS_STD<1e-10)=nan;
        PEAKS_INT=round(PEAKS_MEAN);
        PEAKS_Z=abs(PEAKS_MEAN-PEAKS_INT)./(PEAKS_STD+eps);
        PEAKS_D=abs(PEAKS_MEAN-PEAKS_INT);


        if IS_PLOT
            if III==0
                fprintf('plotting summary graphs [%s %s] ...\n',EXPR,RANGE_msg)
                figure(100+EE);subplot(NSUB1,NSUB2,NSUB2*NSUB1)
                boundedline(0:(NK-1),mean(peaks_num_i),std(peaks_num_i),'k-','transparency', 0.2);hold on;
                plot(0:(NK-1),mean(peaks_num_i),'k-','LineWidth',2);hold on;
                xlabel('Iteration');
                ylabel('Number of peaks');
                title(sprintf ("%s %s",EXPR,RANGE_msg));
                set(gca,'FontSize',14)
            end
            figure(100+EE);subplot(NSUB1,NSUB2,III+1)
            plot_summary_graph(EXPR,xx,HA,INT_RANGE,PEAKS_MEAN,PEAKS_STD,PEAK_SIG,PEAKS_Z,PEAKS_D,PEAK_SIG_TH,LOCSD,PKSD,LOCSA,PKSA);
            title(sprintf ("%s %s",EXPR,RANGE_msg));

        end
        %%
        fprintf('saving final data [%s %s] ...\n',EXPR,RANGE_msg)

        DAT.data=data;
        DAT.raw_data=raw_data;
        DAT.UN=UN;
        DAT.NK=NK;
        DAT.NI=NI;
        DAT.unets=unets;
        DAT.BOOTS=BOOTS;
        DAT.BOOTSi=BOOTSi; %iteration data
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
        %DAT.peaks.peaks_num=peaks_num;
        DAT.peaks.peaks_num_i=peaks_num_i;
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

            %             odfname=sprintf('%s/_raw_data-%s%s.NBOOT-%d.mat',output_dir,EXPR,RANGE_msg,NBOOT);
            %             fprintf('saving raw data to filename (matlab format): %s\n',obfname)
            %             save(odfname,'DAT');

            tifname=sprintf('%s/number_peaks_iterations-%s-%s.NBOOT-%d.txt',output_dir,EXPR,RANGE_msg,NBOOT);
            fprintf('saving table of iteration number bootrapping to: %s\n',tifname)
            mleg=cell(NK,1);
            for ll=1:length(mleg)
                mleg{ll}=sprintf('iter%d',ll);
            end
            ti=array2table(peaks_num_i,'VariableNames',mleg);
            writetable(ti,tifname);

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
%%
% compute peaks for the melodic consonance section
compute_peaks_iterations_melcon_01

