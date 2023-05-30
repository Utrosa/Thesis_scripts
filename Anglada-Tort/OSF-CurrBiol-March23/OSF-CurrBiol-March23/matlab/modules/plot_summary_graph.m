function plot_summary_graph(EXPR,xx,HA,INT_RANGE,PEAKS_MEAN,PEAKS_STD,PEAK_SIG,PEAKS_Z,PEAKS_D,PEAK_SIG_TH,LOCSD,PKSD,LOCSA,PKSA)


%figure(500+EE);clf;
        plot(xx,HA,'-k', 'LineWidth',2);hold on;

        title(EXPR)


        for ll=1:length(PEAKS_MEAN)


            if PEAK_SIG(ll)>PEAK_SIG_TH(1)
                mclr=[1 0 1 0.4];
            elseif PEAK_SIG(ll)>PEAK_SIG_TH(2)
                mclr=[1 0 1 0.2];
            elseif PEAK_SIG(ll)>PEAK_SIG_TH(3)
                mclr=[1 1 0 0.4];
            else
                mclr=[1 1 0 0.2];
            end

            if isnan(PEAKS_MEAN(ll)+PEAKS_STD(ll))
                continue
            end
            rectangle('Position',[PEAKS_MEAN(ll)-PEAKS_STD(ll), 0,PEAKS_STD(ll)*2,PKSD(ll)],'FaceColor',mclr,'EdgeColor',mclr);

        end
        plot(xx,HA,'-k', 'LineWidth',2);hold on;
        plot(LOCSD,PKSD,'cs','MarkerFaceColor','c','MarkerSize',12);hold on;
        plot(LOCSA,PKSA,'bo','MarkerFaceColor','b');hold on;
        for ll=1:length(PEAKS_MEAN)


            text(PEAKS_MEAN(ll),PKSD(ll),sprintf('   %3.2g Z %3.2g D %3.2g',PEAK_SIG(ll),PEAKS_Z(ll),PEAKS_D(ll)));

        end

        for ii=round(min(xx)):round(max(xx))
            plot([ii,ii],[0,max(HA)*1.1],'k--')

        end
        set(gca,'FontSize',14)
        ylim([0,max(HA)*1.1]);
        xlim([min(INT_RANGE)+0.5,max(INT_RANGE)-0.5]);
        set(gca,'XTick',INT_RANGE);