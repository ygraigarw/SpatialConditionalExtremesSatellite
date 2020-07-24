function Prm=PrmS2A(C,Q);
%% function P=PrmS2A(C,Q)
%
% Conditional spatial extremes with delta Laplace residuals
% Philip Jonathan, Rob Shooter
%
% Convert alpha, beta, mu, sigma, delta, rho from structure to array
%
% Input
% *** TO CHECK ***
% C ***
% Q ***
%
% Output
% Prm ? x 1 array

nH=C.nH;
nA=C.nA;

%Q.A is a nH*nA x 1 matrix for nH distances and nA angles for parameters
%- all angles for a given distance first, then all angles for next distance
%Q.A is a 1 x 2 array for limits
%Q.B same as Q.A
%All the others are nH x 1 for parameters and 1 x 2 for limits
if size(Q.A,1)>1; %must be a column vector
   Prm=[reshape(Q.A',nH*nA,1);reshape(Q.B',nH*nA,1);Q.M;Q.S;Q.D;Q.R];
else;
   Prm=[ones(nH*nA,1)*Q.A;ones(nH*nA,1)*Q.B;ones(nH,1)*Q.M;ones(nH,1)*Q.S;ones(nH,1)*Q.D;ones(3,1)*Q.R];
end;
   
return;