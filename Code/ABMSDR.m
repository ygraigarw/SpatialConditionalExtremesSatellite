function [A,B,M,S,D,R,H,BndA,tA,tB]=ABMSDR(X,P,C);
%% function [A,B,M,S,D,R,H,BndA,tA,tB]=ABMSDR(X,P,C);
%
% Conditional spatial extremes with delta-Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross
%
% Calculate alpha, beta, mu, sigma, delta and rho for parameters P and distances H

q=X.q;
H.H0=nan(q,1);
H.A0=nan(q,1);
H.B0=nan(q,1);
for j=1:q;
    %% Use local Cartesian coordinates
    tAvrLtt=(X.Rmt(j,2)+X.Rfr(2))/2;
    tDst=[cos(tAvrLtt*pi/180)*(X.Rmt(j,1)-X.Rfr(1)),X.Rmt(j,2)-X.Rfr(2)];
    H.H0(j)=sqrt(tDst*tDst');
    if tDst(1)>0 && tDst(2)>=0;
       H.A0(j)=atan(tDst(2)/tDst(1));
       H.B0(j)=atan(tDst(2)/tDst(1));
    elseif tDst(1)<0 && tDst(2)>=0;
       H.A0(j)=pi-atan(tDst(2)/(-tDst(1)));
       H.B0(j)=pi-atan(tDst(2)/(-tDst(1)));
    elseif tDst(1)<0 && tDst(2)<0;
       H.A0(j)=pi+atan(tDst(2)/tDst(1));
       H.B0(j)=pi+atan(tDst(2)/tDst(1));
    elseif tDst(1)>0 && tDst(2)<0;
       H.A0(j)=2*pi-atan((-tDst(2))/tDst(1));
       H.B0(j)=2*pi-atan((-tDst(2))/tDst(1));
    elseif tDst(1)==0;
       if tDst(2)>=0;
          H.A0(j)=pi/2;
          H.B0(j)=pi/2;
       else;
          H.A0(j)=3*pi/2;
          H.B0(j)=3*pi/2;
       end;
    end;
end;

%% Critical this ordering matches with definition of H 
for j1=1:q;
   for j2=1:q;
       %% Use local Cartesian coordinates
       tAvrLtt=(X.Rmt(j2,2)+X.Rmt(j1,2))/2;
       tDst=[cos(tAvrLtt*pi/180)*(X.Rmt(j2,1)-X.Rmt(j1,1)),X.Rmt(j2,2)-X.Rmt(j1,2)];
       H.HR(j1,j2)=sqrt(tDst*tDst');
   end;
end;

%% Piecewise linear forms for M, S and D
tH=H.H0;
tDlt=C.Dlt;
nH=size(tH,1);
M=nan(nH,1);
S=nan(nH,1);
D=nan(nH,1);

for i=1:nH;
   tLct=floor(tH(i)/tDlt)+1;
   if tLct<=C.nH; %there are sufficient spline centres
      if rem(tH(i)/tDlt,1)==0;
         M(i)=P.M(tLct);
         S(i)=P.S(tLct);
         D(i)=P.D(tLct);
      elseif tLct+1<=C.nH; 
         tInc=1-rem(tH(i)/tDlt,1);
         M(i)=tInc*P.M(tLct)+(1-tInc)*P.M(tLct+1);
         S(i)=tInc*P.S(tLct)+(1-tInc)*P.S(tLct+1);
         D(i)=tInc*P.D(tLct)+(1-tInc)*P.D(tLct+1);
      else;
         M(i)=P.M(C.nH);
         S(i)=P.S(C.nH);
         D(i)=P.D(C.nH);
      end;
   else; %insufficient spline centres - use constant value
      M(i)=P.M(C.nH);
      S(i)=P.S(C.nH);
      D(i)=P.D(C.nH);
   end;
end;

%% A defined radially parametrically, but angularly in piecewise linear fashion
A=nan(nH,1);
nA=C.nA; %number of rays
tA=nan(nH,nA);
for iA=1:nA;
    for i=1:nH;
        %tA(i,iA)=exp(-(H.H0(i)/P.A1(iA)).^P.A2(iA));
        tLct=floor(tH(i)/tDlt)+1;
        if tLct<=C.nH; %there are sufficient spline centres
            if rem(tH(i)/tDlt,1)==0;
                tA(i,iA)=P.A(tLct,iA);
            elseif tLct+1<=C.nH;
                tInc=1-rem(tH(i)/tDlt,1);
                tA(i,iA)=tInc*P.A(tLct,iA)+(1-tInc)*P.A(tLct+1,iA);
            else;
                tA(i,iA)=P.A(C.nH,iA);
            end;
        else; %insufficient spline centres - use constant value
            tA(i,iA)=P.A(C.nH,iA);
        end;
    end;
end;
tA=[tA tA(:,1)]; %add extra column to correspond with BndA
BndA=(0:nA)'*2*pi/nA; %boundaries 
for i=1:nH;
   for iA=1:nA;
      if H.A0(i)>=BndA(iA) && H.A0(i)<BndA(iA+1);
         tInc=(H.A0(i)-BndA(iA))/(2*pi/nA);
         A(i)=tInc*tA(i,iA)+(1-tInc)*tA(i,iA+1);
         break;
      end;
   end;
end;

%% A defined radially parametrically, but angularly in piecewise linear fashion
B=nan(nH,1);
nB=C.nA; %number of rays
tB=nan(nH,nB);
for iA=1:nA;
    for i=1:nH;
        tLct=floor(tH(i)/tDlt)+1;
        if tLct<=C.nH; %there are sufficient spline centres
            if rem(tH(i)/tDlt,1)==0;
                tB(i,iA)=P.B(tLct,iA);
            elseif tLct+1<=C.nH;
                tInc=1-rem(tH(i)/tDlt,1);
                tB(i,iA)=tInc*P.B(tLct,iA)+(1-tInc)*P.B(tLct+1,iA);
            else;
                tB(i,iA)=P.B(C.nH,iA);
            end;
        else; %insufficient spline centres - use constant value
            tB(i,iA)=P.B(C.nH,iA);
        end;
    end;
end;
tB=[tB tB(:,1)]; %add extra column to correspond with BndA
BndB=(0:nB)'*2*pi/nB; %boundaries 
for i=1:nH;
   for iB=1:nB;
      if H.B0(i)>=BndB(iB) && H.B0(i)<BndB(iB+1);
         tInc=(H.B0(i)-BndB(iB))/(2*pi/nB);
         B(i)=tInc*tB(i,iB)+(1-tInc)*tB(i,iB+1);
         break;
      end;
   end;
end;

R=P.R;

return;