function [P,PNms]=PrmA2S(C,Prm);
%% function P=PrmA2S(Prm)
%
% Conditional spatial extremes with delta-Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross
%
% Convert alpha, beta, mu, sigma, delta, rho from array to structure
%
% Input
% *** TO CHECK ***
% Only use for Prm, never use for Lmt
%
% Output
% P          structure

nH=C.nH;
nA=C.nA;

PNms=cell(2*nH*nA+3*nH+3,1);

j=0;
for i1=1:nH
    for i2=1:nA;
        j=j+1;P.A(i1,i2)=Prm(j);
        PNms{j}=sprintf('A(%g,%g)',i1,i2);
    end;
end;
for i1=1:nH
    for i2=1:nA;
        j=j+1;P.B(i1,i2)=Prm(j);
        PNms{j}=sprintf('B(%g,%g)',i1,i2);
    end;
end;
for i1=1:nH
    j=j+1;
    P.M(i1)=Prm(j);
    PNms{j}=sprintf('M(%g)',i1);
end;
for i1=1:nH
    j=j+1;
    P.S(i1)=Prm(j);
    PNms{j}=sprintf('S(%g)',i1);
end;
for i1=1:nH
    j=j+1;
    P.D(i1)=Prm(j);
    PNms{j}=sprintf('D(%g)',i1);
end;
for i1=1:3
    j=j+1;
    P.R(i1)=Prm(j);
    PNms{j}=sprintf('R(%g)',i1);
end;

return;