function Crr = GetRsdCrr(X,P,H);
%% function Crr = GetRsdCrr(X,P,H)
%
% Conditional spatial extremes with delta-Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross
%
% Overview
% Calculates the correlation matrix for the SCE residual process
% Based on a Gaussian process with powered exponential correlation with distance
% SCE uses the corresponding conditional Gaussian correlation 
%
% Input
% H.H0      q x 1    Distances of conditioned locations from the reference location
% H.HR      q x q    Distances between pairs of conditioned locations
% P.R(1)    1 x 1    Scale parameter for exponential correlation 
% P.R(2)    1 x 1    Exponent parameter for exponential correlation
% X.XRC     n x q    Data for conditioned locations
%
% Output
% Crr       q x q    Correlations for the SCE residual process

q=size(X.XRC,2);

%% START USER INPUT
% Hard-coded parameter suitable for satellite application
% We should make these user input. Just make them fields of ***
Scl1=100; 
Scl2=2;
%END USER INPUT

%% Standard Gaussian correlation "marginally" between conditioned and reference location
Rho=nan(q,1);
for iC=1:q;
    tH=H.H0(iC);
    Rho(iC)=exp( -( tH / (Scl1*(P.R(1))) )^(Scl2*P.R(2)) );
end;

%% Standard Gaussian correlation between pairs of conditioned locations
Crs=nan(q,q);
for iC=1:q;
    for jC=1:q;
        tH=H.HR(iC,jC);
        Crs(iC,jC)=exp( -( tH / (Scl1*(P.R(1))) )^(Scl2*P.R(2)) );
    end;
end;

%% Conditional Gaussian correlation between pairs of conditioned locations
Crr=nan(q,q);
for iC=1:q;
    for jC=1:q;
        t1=sqrt(1-Rho(iC)^2);
        t2=sqrt(1-Rho(jC)^2);
        Crr(iC,jC)=(Crs(iC,jC)-Rho(iC)*Rho(jC))/(t1*t2);
    end;
end;

return;