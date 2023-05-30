addpath('support-functions/'); % dependency on internal libraries
addpath('modules/'); % dependency on modules

close all;clc;clear all;

output_dir='output/';% "To begin at the beginning"/ Under Milk Wood – Dylan Thomas
mkdir(output_dir)

EXPRS1={'btw-2tones'}; % use this data for the singing data
MODEL_DEG=7;


EXPRS2={'melcon-2tones'}; % use this data for preference data

MODELS={'Sing','sim-interval-M','sim-preference-P'}; % the two models
NBOOT=100;

pfname=sprintf('%s/poly_fit_for_manu.%s.txt',output_dir,EXPRS1{1});
prefname=sprintf('%s/pref_smoothes_for_manu.%s.txt',output_dir,EXPRS1{1});

tfname=sprintf('%s/simulations-todo.boot.%d.%s-all.mat',output_dir,NBOOT,EXPRS1{1});


% fixed paramters for all
KernelWidth=0.25; % for 1d kernels
sigma_2d=0.75; % for 2d kernels

XRANGE=[-15,15];
xx=-20:0.1:20;
xxx=min(XRANGE):0.1:max(XRANGE); %% x grid
M=length(xxx);

% for intrval size (motor) model computation.
VAR_TH_RANGE=10; % range to estimate variance where noise it not too noisy

%%
%%%% SINGING DATA %%%%
fprintf('loading singing data!\n')
EE=1;
EXPR=EXPRS1{EE};
[data1,raw_data1,UN1,NK1,NI1,unets1]=analyze_curbio_06(EXPR);

fprintf('Computing marginal last 3 iter... ')
HA=compute_marginal(data1(:,(end-2):end,:),xx,KernelWidth);
HAS=cell(NK1,1);
fprintf('\nComputing marginals... ')
maxx=-1;
for K=1:NK1
    fprintf('.')
    HAS{K}=compute_marginal(data1(:,K,:),xx,KernelWidth);
    maxx=max(maxx,max(HAS{K}));
end
fprintf('\n')

KW=KernelWidth;
xvals=data1(:,1,:);
yvals=data1(:,2,:)-data1(:,1,:); % bias of the first iteration

[pyx_sing,pyv_sing,yvar_sing,yvar2_sing,p1_sing,yp_sing,yvv_sing,yvv2_sing]=compute_biases(xvals,yvals,xx,KW,VAR_TH_RANGE,MODEL_DEG);


fprintf('plotting singing data... ')
figure(10);clf;
for K=1:NK1
    subplot(3,4,K);
    plot(xx,HAS{K},'k','LineWidth',2);hold all;
    title(sprintf('iter=%d',K-1));
    ylim([0 maxx]);
end

figure(11);
clf;subplot(2,1,1);
for I=1:NI1
    plot(xvals(:,1,I),yvals(:,1,I),'.k');hold all;
end

plot(xx,pyx_sing,'k','LineWidth',4);hold all;
plot(xx,xx*0,'g--','LineWidth',2);hold all;

plot(xx,pyv_sing,'b','LineWidth',2);hold all;
plot(xx,yvv_sing,'g:','LineWidth',1)
plot(xx,yvv2_sing,'g:','LineWidth',3)

plot(xx,yp_sing,'r--','LineWidth',3)
xlim(XRANGE);

subplot(2,1,2);
%plot(xx(1:(end-1)),-diff(HA),'k','LineWidth',2)
plot(xx,HA,'k','LineWidth',2);hold on;
xlim(XRANGE);


figure(12);clf;
plot(xx,xx*0,'g--','LineWidth',2);hold all;
xlim(XRANGE);

for I=1:NI1
    plot(xvals(:,1,I),yvals(:,1,I),'.k');hold all;
end
%plot(xvals(:,1,1),yvals(:,1,1),'.k');hold all;
plot(xx,yp_sing,'r-','LineWidth',3);;hold all;

    fprintf('saving temporay files for manu (showing polynomial)')
    pos=abs(xx)<15;
    xxp=xx(pos);
    yp_sing_p=yp_sing(pos);
    writematrix([xxp;yp_sing_p],pfname)
    
%%

%%%% LOAD PREFS DATA AND PLOT IT (COMPUTE SMOOTHING) %%%%

fprintf('loading pref data!\n')

EE=1;
EXPR=EXPRS2{EE};
[data,raw_data,UN,NK,NI,unets]=analyze_curbio_06(EXPR);

fprintf('Computing  smoothed response...')
xvals=data(:,1,1);
yvals=data(:,2,1);
pos=~isnan(xvals+yvals);
xvals=xvals(pos);
yvals=yvals(pos);
pos2=abs(xvals)<max(XRANGE);
xvals=xvals(pos2);
yvals=yvals(pos2);

KW=KernelWidth;
[pyx,pyv]=smooth_dense_data(xvals,yvals,xx,KW);


fprintf('\n')

fprintf('Plotting pref data.\n')


figure(21);clf;
plot(xvals,yvals+0.1*randn(size(yvals)),'.');hold all;
plot(xx,pyx,'k-','LineWidth',2);hold all;
plot(xx,pyv,'c--','LineWidth',2)

figure(22);clf;
plot(xx,pyx,'k-','LineWidth',2);hold all;
for ll=min(XRANGE):max(XRANGE)
    plot([ll,ll],[min(pyx),max(pyx)],'g--');
end

set(gca,'FontSize',13)
ylabel('rating');
xlabel('interval');
set(gca,'XTick',min(XRANGE):max(XRANGE));
xlim(XRANGE);
ylim([min(pyx),max(pyx)]);
title('smoothed prefernce data')

    writematrix([xx;pyx],prefname)


%%

%%%% DO SIMULATIONS %%%%

fprintf('Preparing for simulations:\n')
pyxi=interp1(xx,pyx,xxx); % for pref model in new coordinat system
miu_mapi=interp1(xx,yp_sing,xxx); % for baises in new coordinate system
sig_mapi=interp1(xx,yvv2_sing,xxx);% for baises in new coordinate system

pU=ones(size(xxx,1),1);pU=pU/sum(pU(:)); % uniform distribution

HASi=cell(size(HAS)); % emprical data interpolated to this grid
for K=1:NK1
    HASi{K}=interp1(xx,HAS{K},xxx)';
    HASi{K}=HASi{K}/sum(HASi{K});
end
p0=interp1(xx,HAS{1},xxx)';p0=p0/sum(p0);

% grid search grids
NGRID1=0.1:0.2:4;
NGRID2=0.1:0.2:4;
NGRID3=0.001:0.1:1.201;
NGRID4=0.001:0.1:1.201;


%%% run the two models:
todo=cell(2,1);
todo{1}.NGRID_A=NGRID3;
todo{1}.NGRID_B=NGRID4;
todo{1}.todo='M';
todo{1}.name='Interval-size-model';
todo{1}.var_names={'noise-factor','bias-factor'};
todo{1}.clr='b';

todo{2}.NGRID_A=NGRID1;
todo{2}.NGRID_B=NGRID2;
todo{2}.todo='P';
todo{2}.name='Preference-model';
todo{2}.var_names={'gamma','sigma'};
todo{2}.clr='r';


NT=length(todo);

% prepare all figures:
figure(31);clf;figure(32);clf;figure(41);clf;figure(42);clf;figure(51);clf;figure(52);clf;figure(61);clf;figure(62);clf;



fprintf('Computing simulations:\n')
HA0=compute_marginal(data1(:,(end-2):end,:),xxx,KernelWidth); % singing data marginal
for T=1:length(todo)
    NGRID_A=todo{T}.NGRID_A;
    NGRID_B=todo{T}.NGRID_B;
    figure(30+T);clf;

    Ascores=nan(length(NGRID_A),length(NGRID_B),NBOOT);
    Bscores=nan(length(NGRID_A),length(NGRID_B));
    for I1=1:length(NGRID_A)
        fprintf("\n*")
        for I2=1:length(NGRID_B)
            fprintf('+')
            HA1s=nan(NBOOT,length(HA0));
            HA0s=nan(NBOOT,length(HA0));
            Ascores_boot=nan(1,NBOOT);
            parfor B=1:NBOOT

                if B==1
                    idx=1:size(data1,1);
                else
                    idx=randi(size(data1,1),1,size(data1,1)  );
                end

                % real data
                sing_data=data1(idx,:,:);
                HAS0=compute_marginal(sing_data(:,(end-2):end,:),xx,KernelWidth);
                HAS=cell(NK1,1);
                HASi=cell(size(HAS)); % emprical data interpolated to this grid
                %fprintf('\nComputing marginals... ')
                for K=1:NK1
                    %   fprintf('.')
                    HAS{K}=compute_marginal(sing_data(:,K,:),xx,KernelWidth);
                    HASi{K}=interp1(xx,HAS{K},xxx)';
                    HASi{K}=HASi{K}/sum(HASi{K});
                end
                p1=nan;
                p2=nan;
                switch todo{T}.todo

                    case 'M'
                        %NOISE_COMBINED_REDUCE_FACTOR
                        %BIAS_COMBINED_REDUCE_FACTOR
                        p1=NGRID_A(I1);
                        p2=NGRID_B(I2);
                    case 'P'
                        %gamma
                        %sigma_sen
                        p1=NGRID_A(I1); % params
                        p2=NGRID_B(I2); % params
                    otherwise
                        fprintf('Wrong todo value\n');
                        assert(1==0);
                end
                [sim_data,HAS_sim]=do_models(todo{T}.todo,p1,p2,sing_data,NK1,xxx,miu_mapi,sig_mapi,pyxi,p0,KernelWidth);


                HA1=compute_marginal(sim_data(:,(end-2):end,:),xxx,KernelWidth);
                rscore1=JSD2(HA0,HA1);


                HA1s(B,:)=HA1;
                HA0s(B,:)=HA0;
            end
            Ascores(I1,I2,:)=Ascores_boot;
            Bscores(I1,I2)=JSD2(mean(HA0s,1),mean(HA1s,1));
            figure(50+T);
            imagesc(NGRID_B,NGRID_A,Bscores);colorbar;hold on;axis xy
            [i1,i2]=find(Bscores==min(Bscores(:)),1);
            p2=NGRID_B(i2);
            p1=NGRID_A(i1);

            if (I1==i1) && (I2==i2) % new best
                figure(60+T);clf;
                boundedline(xxx,HA0s(1,:),std(HA0s,1),'k-','transparency', 0.1);hold on;
                boundedline(xxx,HA1s(1,:),std(HA1s,1),'r-','transparency', 0.1);hold on;
                set(gca,'FontSize',14);
                maxf=max(max(HA0s(:),HA1s(:)));
                ylim([0,maxf]);
            end


            plot(NGRID_B(i2),NGRID_A(i1),'ow');hold on;

            msg=sprintf('%s %s=%g %s=%g',todo{T}.name,todo{T}.var_names{1},NGRID_A(i1),todo{T}.var_names{2},NGRID_B(i2));
            title(msg);

            ylabel(todo{T}.var_names{1});
            xlabel(todo{T}.var_names{2});

            drawnow
            %todo{T}.Ascores_red=Ascores_red;
            todo{T}.Ascores=Ascores;
            todo{T}.Bscores=Bscores;
            todo{T}.Best=[NGRID_A(i1),NGRID_B(i2)];


        end
    end

end


fprintf("\n")
fprintf('saving to target file all workspace %s\n',tfname)
save(tfname);

%%
tic()
%%%% BOOTRAP TH FINAL VALUE AND PLOT THINGS %%%%
fprintf('Computing simulations bootrapping:\n')
HA0=compute_marginal(data1(:,(end-2):end,:),xxx,KernelWidth); % singing data marginal

for T=1:length(todo)
    NGRID_A=todo{T}.NGRID_A;
    NGRID_B=todo{T}.NGRID_B;
    p1=todo{T}.Best(1);
    p2=todo{T}.Best(2);

    HA1s=nan(NBOOT,length(HA0));
    HA0s=nan(NBOOT,length(HA0));
    for B=1:NBOOT
        if mod(B,10)==1
            fprintf('+')
        end

        if B==1
            idx=1:size(data1,1);
        else
            idx=randi(size(data1,1),1,size(data1,1)  );
        end

        % real data
        sing_data=data1(idx,:,:);
        HAS0=compute_marginal(sing_data(:,(end-2):end,:),xx,KernelWidth);
        HAS=cell(NK1,1);
        HASi=cell(size(HAS)); % emprical data interpolated to this grid

        for K=1:NK1
            HAS{K}=compute_marginal(sing_data(:,K,:),xx,KernelWidth);
            HASi{K}=interp1(xx,HAS{K},xxx)';
            HASi{K}=HASi{K}/sum(HASi{K});
        end


        [sim_data,HAS_sim]=do_models(todo{T}.todo,p1,p2,sing_data,NK1,xxx,miu_mapi,sig_mapi,pyxi,p0,KernelWidth);


        HA1=compute_marginal(sim_data(:,(end-2):end,:),xxx,KernelWidth);
        rscore1=JSD2(HA0,HA1);
        Ascores(I1,I2,B)=rscore1;
        HA1s(B,:)=HA1;
        HA0s(B,:)=HA0;

        if B==1
            figure(40+T);clf;
            plot_simulations(HASi,HAS_sim,xxx,XRANGE,todo{T}.name);
            HASi_tosave=nan(size(HASi,1),length(HASi{1}));
            for K=1:NK
                HASi_tosave(K,:)=HASi{K};
            end

            HAS_sim_tosave=nan(size(HAS_sim,1),length(HAS_sim{1}));
            for K=1:NK
                HAS_sim_tosave(K,:)=HAS_sim{K};
            end

            fname1=sprintf('%s/simulations-sample.%d.%s.csv',output_dir,NBOOT,EXPRS1{1});
            fname2=sprintf('%s/simulations-sample.%d.%s-%s.csv',output_dir,NBOOT,EXPRS1{1},todo{T}.name);
            writematrix([xxx;HASi_tosave],fname1)
            writematrix([xxx;HAS_sim_tosave],fname2)

        end

    end
    fprintf('\n')
    Bscores(I1,I2)=JSD2(mean(HA0s,1),mean(HA1s,1));
    figure(60+T);clf;
    boundedline(xxx,HA0s(1,:),std(HA0s,1),'k-','transparency', 0.1);hold on;
    boundedline(xxx,mean(HA1s(:,:)),std(HA1s,1),'r-','transparency', 0.1);hold on;
    set(gca,'FontSize',14);
    maxf=max(max(HA0s(:),HA1s(:)));
    ylim([0,maxf]);

    figure(65+T);clf;
    boundedline(xxx,HA0s(1,:),std(HA0s,1),'k-','transparency', 0.1);hold on;
    boundedline(xxx,HA1s(1,:),std(HA1s,1),'r-','transparency', 0.1);hold on;
    set(gca,'FontSize',14);
    maxf=max(max(HA0s(:),HA1s(:)));
    ylim([0,maxf]);


    todo{T}.BEST.HA0s=HA0s;
    todo{T}.BEST.HA1s=HA1s;

    figure(70+T);clf;
    idx=1:size(data1,1);
    sing_data=data1(idx,:,:);
    HAS0=compute_marginal(sing_data(:,(end-2):end,:),xx,KernelWidth);
    HAS=cell(NK1,1);
    HASi=cell(size(HAS)); % emprical data interpolated to this grid
    for K=1:NK1
        HAS{K}=compute_marginal(sing_data(:,K,:),xx,KernelWidth);
        HASi{K}=interp1(xx,HAS{K},xxx)';
        HASi{K}=HASi{K}/sum(HASi{K});
    end

    [sim_data,HAS_sim]=do_models(todo{T}.todo,p1,p2,sing_data,NK1,xxx,miu_mapi,sig_mapi,pyxi,p0,KernelWidth);
    HA1=compute_marginal(sim_data(:,(end-2):end,:),xxx,KernelWidth);
    plot_simulations(HASi,HAS_sim,xxx,XRANGE,todo{T}.name)


    fname1=sprintf('%s/simulations-runs.boot.%d.%s.csv',output_dir,NBOOT,EXPRS1{1});
    fname2=sprintf('%s/simulations-runs.boot.%d.%s-%s.csv',output_dir,NBOOT,EXPRS1{1},todo{T}.name);

    fprintf('Saving in time %s',fname1);
    writematrix([xxx;HA0s],fname1);

    fprintf('Saving in time %s',fname2);
    writematrix([xxx;HA1s],fname2);



    fprintf('\nComputing kernels\n')
    fprintf('KERNELS:real data:      :\t')
    [FA,FAS]=compute_kernels(data1,xxx,xxx,sigma_2d);

    fprintf('KERNELS:simualted data:\t')
    [FA_sim,FAS_sim]=compute_kernels(sim_data,xxx,xxx,sigma_2d);



    fprintf('plotting all kernels... \n')

    figure(70+T);clf
    NSUB1=4;
    NSUB2=6;

    myKS=1:2:NK1;

    MX=max(max(FAS{K}(2:end-1,2:end-1)));

    for KK=1:length(myKS)

        K=myKS(KK);
        subplot(NSUB1,NSUB2,KK)
        imagesc(xxx(2:end-1),xxx(2:end-1),FAS{K}(2:end-1,2:end-1));axis xy;
        title(sprintf('singing data iter=%d',K-1),'Color','k')

        subplot(NSUB1,NSUB2,KK+1*NSUB2)
        imagesc(xxx(2:end-1),xxx(2:end-1),FAS_sim{K}(2:end-1,2:end-1));axis xy;
        title(sprintf('sim %s iter=%d',todo{T}.name,K-1),'Color',todo{T}.clr)

    end
end

toc()
tfname2=sprintf('%s/temp-final-simulations-todo.boot.%d.%s-all.mat',output_dir,NBOOT,EXPRS1{1});
fprintf('temporary last --> saving to target file all workspace %s\n',tfname)
save(tfname2)
