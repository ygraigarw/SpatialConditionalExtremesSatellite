function PltDgn(C,iI);
%% function PltDgn(C,iI);
%
% Conditional spatial extremes with delta-Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross
%
% Conditional quantile constraints (Keef et al., 2013) diagnostic plots,
% and MCMC diagnostics (acceptance rate and proposal standard deviation)

clf;
q=size(C.Alp,2);
tClr=jet(q);

subplot(1,2,1);hold on;
for iL=1:q
    plot(C.Alp(max(1,iI-9):iI,iL),C.Bet(max(1,iI-9):iI,iL),'--','color',tClr(iL,:));
    plot(C.Alp(iI,iL),C.Bet(iI,iL),'.','color',tClr(iL,:),'markersize',20); title('\beta on \alpha');
    if C.IsCQC==1
        plot(C.CQCBndL{iL}(:,1),C.CQCBndL{iL}(:,2),'color',tClr(iL,:));
    end
end
HlpFnt; HlpAxsLmt;
set(gca,'xlim',[min(C.Alp(iI,:))-0.05,max(C.Alp(iI,:))+0.05]);
set(gca,'ylim',[min(C.Bet(iI,:))-0.05,max(C.Bet(iI,:))+0.05]);
title('\beta on \alpha');

subplot(2,2,2); 
plot(C.Ngt(1:iI,:)); title('Proposal s.d.'); 
HlpFnt; HlpAxsLmt;

subplot(2,2,4); 
plot(C.AccRat(1:iI,:)); title('Acceptance rate'); 
HlpFnt; HlpAxsLmt;

drawnow;

return;