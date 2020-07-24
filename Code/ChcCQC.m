function isOK=ChcCQC(lagAlp,lagBet,lagMu,psi,X);
%% function isOK=ChcCQC(lagAlp,lagBet,lagMu,psi,X);
%
% Conditional spatial extremes with delta-Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross
%
% Check that conditional quantile constraints (Keef et al., 2013), if used,
% are not violated for the given parameters

nLag=size(lagAlp,1);

for iL=1:nLag
   
   x=X.X0C;
   y=X.XRC(:,iL);
   
   isOK=CQC(lagAlp(iL),lagBet(iL),lagMu(iL),psi,x,y);
   if isOK==0
      break;
   end
   
end;

return;