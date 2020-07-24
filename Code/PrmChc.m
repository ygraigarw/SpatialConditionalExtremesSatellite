function IsOK=PrmChc(Prm,PrmLmt);
%% function IsOK=PrmChc(Prm,PrmLmt);
%
% Conditional spatial extremes with delta-Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross
%
% Check that parameter values are valid (i.e., within limits)
% All inputs as arrays, not structures

IsOK=0;

t1=Prm<PrmLmt(:,1);
t2=Prm>PrmLmt(:,2);

if sum(t1+t2)==0;
   IsOK=1;
end;

return;