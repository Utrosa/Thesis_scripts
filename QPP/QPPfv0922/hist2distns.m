function lgd=hist2distns(SamplesRef, Samples,SampleRefName, SampleName, fc_ref,fc_sample, bin)

h2null=histogram(SamplesRef,'Normalization','probability', 'BinWidth', bin,'FaceColor',fc_ref); hold on; %SystematicStim Null
h2=histogram(Samples,'Normalization','probability', 'BinWidth', bin,'FaceColor',fc_sample,'FaceAlpha',.4); %axis([0 max([QPPdelay3All; QPPdelay3restAll]) 0 1])
h2BinCenters = h2.BinEdges + h2.BinWidth/2; h2BinCenters=h2BinCenters(1:length(h2.Values));
h2nullBinCenters = h2null.BinEdges + h2null.BinWidth/2; h2nullBinCenters=h2nullBinCenters(1:length(h2null.Values));
h2Values=h2.Values; h2nullValues=h2null.Values;
if length(h2.Values)<length(h2null.Values)    
    h2Values(length(h2Values)+1:length(h2nullValues))=0;
    h2BinCenters=h2nullBinCenters;
elseif length(h2.Values)>length(h2null.Values)
    h2nullValues(length(h2nullValues)+1:length(h2Values))=0;
else
    h2BinCenters=h2nullBinCenters;
end
lgd=legend({[SampleRefName ' (m=' num2str(mean(SamplesRef)) ', ' num2str(length(SamplesRef)) 'pts)'], [SampleName ' (m=' num2str(mean(Samples)) ', ' num2str(length(Samples)) 'pts)']});
lgd.FontSize = 8;
% kl = kldiv(h2BinCenters,h2nullValues+eps, h2Values+eps);
% js = kldiv(h2BinCenters,h2nullValues+eps, h2Values+eps,'js');
% ks = kstest2(SamplesRef,Samples,'Alpha',ksalpha); 