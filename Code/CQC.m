function isOK=CQC(alp,bet,mu,psi,x,y);
%% function isOK=CQC(alp,bet,mu,psi,x,y);
%
% Conditional spatial extremes with delta-Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross
%
% Calculates conditional quantile constraint boundaries (Keef et al., 2013)
% and confirms compliance

% INPUT *** TO CHECK ***
% x : threshold exceedances for reference location
% y : corresponding values at remote location

% OUTPUT
% isOK : =0, violates constraints
% isOK : =1, satisfies constraints

isOK=0;

%Conditional quantile constraints
z=max((y-alp*x)./x.^bet);
zplus=max(y-x);
nu=max(x);

C1=1;
C2=1-bet*z*nu^(bet-1);
C3=1-(nu^(bet-1))*z + zplus/nu;
T1=alp<=min([C1;C2;C3]);

C4=1-bet*z*nu^(bet-1);
T2=alp>C4;

C5=(1-1/bet)*((bet*z)^(1/(1-bet)))*((1-alp).^(-bet/(1-bet)))+zplus;
T3=C5>0;

if isempty(mu*psi)==0
   %additional constraint of Lugrin (2018) - but not needed in practice
   C6=alp+bet*(psi^(bet-1))*mu;
   T4=C6>0;
else
   T4=1; %constraint always satisfied
end

if T4>0 %Lugrin (2018) constraint MUST be satisfied
   if T1>0 %We could have T1>0 OR T2*T3>0
      isOK=1;
   elseif T2*T3>0
      isOK=1;
   end
end

return;