function Nll = SCENll(P,C,X)
%% function Nll = SCENll(P,C,X)
%
% Conditional spatial extremes with delta-Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross
%
% Input
% *** TO CHECK ***
%
% Output
% Nll         scalar negative log likelihood from delta Laplace

%% Evaluate parameters at obervations
[A,B,M,S,D,~,H]=ABMSDR(X,P,C); %all parameters are q x 1

%% Define threshold based on non-exceedance probability
psi=X.ThrLpl;

%% Optionally check Keef constraints
if C.IsCQC==1
   isOK=ChcCQC(A,B,M,psi,X);
   if isOK==0
      Nll=inf;
      return;
   end
end

%% Set-up
[n,q]=size(X.XRC); %number of observations and remote locations
PdfZ=nan(n,q);CdfZ=nan(n,q);QntNrm=nan(n,q);PdfNrm=nan(n,q); %initialise arrays

%% Estimate the covariance matrix (\Sigma in text) and related quantities
Crr=GetRsdCrr(P,H,C);
invCrr=inv(Crr);
invCrrSR=sqrtm(invCrr);
lgrDtrCrr=log(abs(det(Crr)));
if isreal(invCrrSR)==0
    Nll=inf;
    return;
end

%% Evaluate nll
Mth=1; %=1 == delta-Laplace; =2 == Gaussian (for checking only)
switch Mth

    case 1 %delta-Laplace
        
        % pdf and cdf
        for j=1:q
            Kpp=sqrt(gamma(1/D(j))/gamma(3/D(j)));
            tM=X.X0C*A'+(X.X0C*ones(1,q)).^(ones(n,1)*B').*(ones(n,1)*M');
            tS=(X.X0C*ones(1,q)).^(ones(n,1)*B').*(ones(n,1)*S')*Kpp;
            PdfZ(:,j)=(D(j)./(2*tS(:,j).*gamma(1/D(j)))).*exp(-abs((X.XRC(:,j)-tM(:,j))./tS(:,j)).^D(j));
            CdfZ(:,j)=0.5+0.5*sign(X.XRC(:,j)-tM(:,j)).*gammainc(abs((X.XRC(:,j)-tM(:,j))./tS(:,j)).^D(j),1/D(j));
            QntNrm(:,j)=norminv(CdfZ(:,j));
            PdfNrm(:,j)=normpdf(QntNrm(:,j));
        end %gammainc as defined by MATLAB - careful!
        
        % evaluate the negative log likelihood
        tnll=nan(n,1);
        for i=1:n
            tnll(i)=(q/2)*log(2*pi)+(1/2)*lgrDtrCrr+(1/2)*sum((invCrrSR*QntNrm(i,:)').^2)-sum(log(PdfZ(i,:)./PdfNrm(i,:))); %need to check this is faster
        end
        Nll=sum(tnll);
        if isnan(Nll)==1
            Nll=inf;
        end
        
    case 2 %Gaussian (useful check)
        
        % evaluate the negative log likelihood
        tnll=nan(n,1);
        for i=1:n
            tQntNrm=(Z(i,:)'-M)./S;
            tnll(i)=(q/2)*log(2*pi)+sum(log(S))+(1/2)*lgrDtrCrr+(1/2)*sum((invCrrSR*tQntNrm).^2);
        end
        Nll=sum(tnll);
        if isnan(Nll)==1
            Nll=inf;
        end
        
end

return;