function PltPrm(C,nI,nToPlt,tStr)
%% function PltPrm(C,nI,nToPlt,tStr);
%
% Conditional spatial extremes with delta-Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross
%
% Plots of SCE model parameters with direction and distance

%Handle whether function is generating intermediate or final plots
if nargin==3
    tStr=[];
end

nDlt=min(nI,nToPlt);

Prm=C.Prm(nI-nDlt+1:nI,:);
medP=PrmA2S(C,median(Prm)');
lowP=PrmA2S(C,quantile(Prm,0.025)');
uppP=PrmA2S(C,quantile(Prm,0.975)');

Dst=(0:C.Dlt:C.HMxm)';

figure(3); clf;

subplot(2,3,1); hold on;
if C.nA==1
    plot(Dst,lowP.A,'k--');
    plot(Dst,uppP.A,'k--');
    plot(Dst,medP.A,'k-');
else
    for iA=1:C.nA
        plot(Dst,medP.A(:,iA),'color',pClr(iA));
    end
end
HlpFnt; HlpAxsLmt;
title '\alpha';
subplot(2,3,2); hold on;
if C.nA==1
    plot(Dst,lowP.B,'k--');
    plot(Dst,uppP.B,'k--');
    plot(Dst,medP.B,'k-');
else
    for iA=1:C.nA
        plot(Dst,medP.B(:,iA),'color',pClr(iA));
    end
end
HlpFnt; HlpAxsLmt;
title '\beta';
subplot(2,3,3); hold on;
plot(Dst,lowP.M,'k--');
plot(Dst,uppP.M,'k--');
plot(Dst,medP.M,'k-');
HlpFnt; HlpAxsLmt;
title '\mu';
subplot(2,3,4); hold on;
plot(Dst,lowP.S,'k--');
plot(Dst,uppP.S,'k--');
plot(Dst,medP.S,'k-');
HlpFnt; HlpAxsLmt;
title '\sigma';
subplot(2,3,5); hold on;
plot(Dst,lowP.D,'k--');
plot(Dst,uppP.D,'k--');
plot(Dst,medP.D,'k-');
HlpFnt; HlpAxsLmt;
title '\delta';
subplot(2,3,6); 
histogram(Prm(:,end-2), 'Normalization','pdf', 'FaceColor','blue')
hold on
histogram(Prm(:,end-3), 'Normalization','pdf', 'FaceColor','red')
title '\rho';
legend({'\rho_1','\rho_2'})
%title 
HlpFnt; HlpAxsLmt;
title '\rho';
if isempty(tStr)==0 
    HlpSveImg(tStr{1},2); 
end

% Set number of rows and columns per plot
RC=[1 1;
    1 2;
    1 3;
    2 2;
    2 3;
    2 3;
    2 4;
    2 4;
    3 3;
    3 4;
    3 4;
    3 4;
    3 5;
    3 5;
    3 5];


PltLim = [ min([lowP.A(:);lowP.B(:)]), max([uppP.A(:);uppP.B(:)]) ];

figure(4); clf;
for iA=1:C.nA
    subplot(RC(C.nA,1),RC(C.nA,2),iA); hold on;
    plot(Dst,lowP.A(:,iA),'k--');
    plot(Dst,uppP.A(:,iA),'k--');
    plot(Dst,medP.A(:,iA),'k-');
    HlpFnt; HlpAxsLmt; set(gca,'ylim',PltLim);
    title(sprintf('\\alpha(:,%g)',iA));
end
if isempty(tStr)==0 
    HlpSveImg(tStr{2},2); 
end

figure(5); clf;
for iA=1:C.nA
    subplot(RC(C.nA,1),RC(C.nA,2),iA); hold on;
    plot(Dst,lowP.B(:,iA),'k--');
    plot(Dst,uppP.B(:,iA),'k--');
    plot(Dst,medP.B(:,iA),'k-');
    HlpFnt; HlpAxsLmt; set(gca,'ylim',PltLim);
    title(sprintf('\\beta(:,%g)',iA));
end
if isempty(tStr)==0
    HlpSveImg(tStr{3},2); 
end

return;