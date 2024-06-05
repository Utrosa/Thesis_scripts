function plotNets(ROI2Net,NetLB, PL,net)
if isempty(PL), PL=50; end

PLhs=round(PL/2)+1; PLhe=round(PL/2)+PL;
nnet=length(unique(ROI2Net)); iROI2NetC=cell(length(NetLB),1); iNetL=zeros(nnet+1,1);
[a,iROI2Net]=sort(ROI2Net); NetID=unique(ROI2Net)';
for inet=NetID
    iROI2NetC{inet}=iROI2Net(a==inet);
    iNetL(inet+1)=iNetL(inet)+length(iROI2NetC{inet});
end
nX2=iNetL(end); tick_pos_msk=diff(iNetL)/2+iNetL(1:end-1);
hold on
for i=2:length(NetLB)        
    plot([0 nX2*2],[iNetL(i) iNetL(i)],'k'); 
end             
set(gca,'XTick',[PLhs (PLhs+PLhe)/2 PLhe],'XTickLabel',[1 round(PL/2) PL],'YTick',tick_pos_msk,'YTickLabel',NetLB)

if net==1    
    for i=2:length(NetLB)        
        plot([iNetL(i) iNetL(i)],[0 nX2],'k'); 
    end        
    axis equal; axis([0 nX2 0 nX2]); 
    set(gca,'XTick',tick_pos_msk,'XTickLabel',NetLB)
    xtickangle(45)
end
hold off;