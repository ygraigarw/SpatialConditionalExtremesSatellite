function PltTrc(C,iI,nToPlt);
%% function PltTrc(C,iI,nToPlt);
%
% Conditional spatial extremes with delta Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross
%
% MCMC trace plots

nPlt=size(C.Prm,2)+3;
[~,PNms]=PrmA2S(C,C.Prm(iI,:)');

strI=max(0,iI-nToPlt)+1;

if rem(sqrt(nPlt),1)==0
   td1=sqrt(nPlt);
   td2=sqrt(nPlt);
else
   td1=floor(sqrt(nPlt));
   td2=floor(sqrt(nPlt))+2;
end;
for j=1:size(C.Prm,2);
   subplot(td1,td2,j); hold on;
   plot(C.Prm(strI:iI,j)); title(PNms{j}); HlpFnt; HlpAxsLmt;
end;
subplot(td1,td2,nPlt-2); plot((1:iI)',C.Nll(1:iI)); title('NLL'); HlpFnt; HlpAxsLmt;
subplot(td1,td2,nPlt-1); plot((strI:iI)',C.Nll(strI:iI)); title('NLL'); HlpFnt; HlpAxsLmt;
subplot(td1,td2,nPlt); plot((iI-9:iI)',C.Nll(iI-9:iI)); title('NLL'); HlpFnt; HlpAxsLmt;
drawnow;

return;