function C=MCMC(X,C)
%% function C=MCMC(X,C);
%
% Conditional spatial extremes with delta-Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross

figure(2); clf;   

%% Specify parameter limits 
PrmLmt=PrmS2A(C,C.Lmt); %parameter limits as numeric array (converting from string)
nP=size(PrmLmt,1); %number of parameters

%% Local copies of parameters (for brevity)
q=X.q; %number of remote locations
nI=C.nI; %number of iterations
AdpItr1=C.AdpItr1; %number of iterations for fixed nugget
AdpItr2=C.AdpItr2; %additional iterations for adaptive nugget
AdpBet=C.AdpBet; %adaptive MCMC parameter from Roberts and Rosenthal (2009)
NmbToPlt=C.NmbToPlt; %number of iterations from end to plot
NgtStr=C.NgtStr; %starting nuggest standard deviation

%% New analysis, delete any previous progress files (BUT KEEP CQC BOUNDARIES)
if C.StrAgn==1
    fprintf('Starting again. Deleting previous analysis.\n');
    tFil=sprintf('%s-Prm%.2f.mat',X.DatNam,100*X.Nep);
    if exist(fullfile(cd,tFil),'file')==2
        tStr=sprintf('delete %s\n',tFil);
        eval(tStr);
    end
end

%% Handle restarts
tFil=sprintf('%s-Prm.mat',X.DatNam);
if exist(fullfile(cd,tFil),'file')==2
    fprintf(1,'Using previous chain\n');
    load(tFil,'C');
    tC=C; 
    % make sure mI is correct, when restarting from deliberately crashed runs
    mI=size(tC.Alp,1);
    tC.AccRat=tC.AccRat(1:mI,:); % previous acceptance rate values
    tC.Ngt=tC.Ngt(1:mI,:); %previous nugget values
    tC.Nll=tC.Nll(1:mI,:); %previous negative log-likelihood values
    tC.Prm=tC.Prm(1:mI,:); %previous SCE model parameter values
    C = rmfield(C,{'AccRat','Ngt','Nll','Prm','nI'});

    % MCMC arrays
    C.AccRat=nan(mI+nI,nP);C.AccRat(1:mI,:)=tC.AccRat; %acceptance rate
    C.Ngt=nan(mI+nI,nP);C.Ngt(1:mI,:)=tC.Ngt; %adaptive MCMC nugget
    C.Nll=nan(mI+nI,1);C.Nll(1:mI,:)=tC.Nll; %negative log-likelihood
    C.Prm=nan(mI+nI,nP);C.Prm(1:mI,:)=tC.Prm; %SCE model parameter values
    C.nI=mI+nI; %number of iterations

    Nll=C.Nll(mI);
    Prm=C.Prm(mI,:)';
else
    mI=0;
    %MCMC arrays
    C.AccRat=nan(nI,nP);C.AccRat(1,:)=0; %acceptance rate
    C.Ngt=nan(nI,nP); %adaptive MCMC nugget
    C.Nll=nan(nI,1); %negative log-likelihood
    C.Prm=nan(nI,nP); %SCE model parameter values
end

%% Loop over MCMC iterations
%
fprintf(1,'Starting MCMC analysis for case NEP = %g.\n',X.Nep); %print progress to terminal 
tic; %timer (returns time for each 100 iterations)

for iI=mI+1:mI+nI  %loop over iterations
    
    %find valid starting solution if iI=1
    if iI==1
        PrmStr=C.Prm0;
        tP=PrmA2S(C,PrmStr);
        tP.p=q+1;
        NllStr=SCENll(tP,C,X); %calculate NLL of proposed solution
        %
        if isinf(NllStr)==1 %no valid starting solution found. Terminate.
            fprintf(1,'Warning: invalid starting solution. Terminating.\n');
            return;
        else %make the current state the starting state for MCMC
            Prm=PrmStr;
            Nll=NllStr;
            fprintf(1,'Starting solution found. Starting MCMC\n');
        end
    end
    
    % iteration counter on screen
    if rem(iI,10)==0
        fprintf(1,'+');
    else
        fprintf(1,'.');
    end
    if rem(iI,100)==0
        fprintf(1,'\n');
    end
    
    % loop over parametric forms
    for iP=1:nP %Metropolis Hastings in Gibbs, one parameter at a time
        
        % define candidate in terms of current state
        PrmC=Prm;
        
        nH=C.nH; nA=C.nA;
        if iI<=AdpItr1 %fixed nugget
            tNgt=NgtStr;
            PrmC(iP)=PrmC(iP)+randn*tNgt;
        elseif iI<=AdpItr2 %adjust nugget standard devation
            if rem(iI,1)==0
                if C.AccRat(iI-1,iP)>0.3
                    tNgt=C.Ngt(iI-1,iP)*1.01;
                elseif C.AccRat(iI-1,iP)<0.2
                    tNgt=C.Ngt(iI-1,iP)*0.99;
                else
                    tNgt=C.Ngt(iI-1,iP);
                end
                if tNgt>2 %max proposal s.d. is hard coded to TWO
                    tNgt=2;
                end
            else
                tNgt=C.Ngt(iI-1,iP);
            end
            PrmC(iP)=PrmC(iP)+randn*tNgt;
        else %adaptive Metropolis for A, B, M, S
            if iP<=nH %update alpha, beta, mu and sigma together
                jP=[nHnA*(iP-1)+(1:nHnA), nHnA*(nH+iP-1)+(1:nHnA), 2*nHnA*nH+iP, 2*nHnA*nH+nH+iP]';   % Following lines are breakdown of equation XX in Roberts & Rosenthal (2009)
                nJ=size(jP,1);
                t1=real((1-AdpBet)*2.38*sqrtm(cov(C.Prm(max(1,iI-999):iI-1,jP))/nJ)*randn(nJ,1));  
                t2=AdpBet*0.1*(randn(nJ,1)/nJ);
                PrmC(jP)=PrmC(jP)+t1+t2;
                tNgt=NaN;
                C.Ngt(iI-1,jP(2:end))=NaN;
                C.AccRat(iI-1,jP(2:end))=NaN;
            elseif iP>2*nHnA*nH+2*nH
                if rem(iI,1)==0
                    if C.AccRat(iI-1,iP)>0.3
                        tNgt=C.Ngt(iI-1,iP)*1.01;
                    elseif C.AccRat(iI-1,iP)<0.2
                        tNgt=C.Ngt(iI-1,iP)*0.99;
                    else
                        tNgt=C.Ngt(iI-1,iP);
                    end
                    if tNgt>2 %max proposal sd is hard coded to TWO
                        tNgt=2;
                    end
                else
                    tNgt=C.Ngt(iI-1,iP);
                end
                PrmC(iP)=PrmC(iP)+randn*tNgt;
            end
        end
        
        % evaluate likelihood at current state and candidate
        PC=PrmA2S(C,PrmC);
        IsOK=PrmChc(PrmC,PrmLmt);
        if IsOK==1
            NllC=SCENll(PC,C,X);
        else
            NllC=inf;
        end
        
        % MH acceptance step
        if (exp(-NllC+Nll) > rand) && isinf(NllC)==0 && isnan(NllC)==0;
            Prm=PrmC;
            Nll=NllC;
            if iI>1
                if iI>100 %only use last 100 iterations to adjust acceptance rate
                    jI=100;
                else
                    jI=iI;
                end
                if iI<=AdpItr2
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1)+1)/jI;
                else
                    if iP<=nH
                        C.AccRat(iI,jP(1))=(C.AccRat(iI-1,jP(1))*(jI-1)+1)/jI;
                    elseif iP>2*nHnA*nH+2*nH
                        C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1)+1)/jI;
                    end
                end
            end
        else
            if iI>1 %only use last 100 iterations to adjust acceptance rate
                if iI>100
                    jI=100;
                else
                    jI=iI;
                end
                if iI<=AdpItr2
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1))/jI;
                else
                    if iP<=nH
                        C.AccRat(iI,jP(1))=(C.AccRat(iI-1,jP(1))*(jI-1))/jI;
                    elseif iP>2*nHnA*nH+2*nH
                        C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1))/jI;
                    end
                end
            end
        end
        
        % save proposal standard deviation
        C.Ngt(iI,iP)=tNgt;
        
    end  %end loop over MCMC iterations
    
    % update after complete iteration over variables
    C.Prm(iI,:)=Prm';
    C.Nll(iI)=Nll;
    [tA,tB]=ABMSDR(X,PrmA2S(C,C.Prm(iI,:)'),C);
    C.Alp(iI,:)=tA';
    C.Bet(iI,:)=tB';
    
    %upon start-up (in the first 50 iterations), update trace plots after every 10th iteration
    if iI<=50 
        if rem(iI,10)==0 
            % obtain trace plots
            figure(1); clf;
            PltTrc(C,iI,iI);
            tFil=sprintf('%s-TrcPrg%.2f',X.DatNam,100*X.Nep);
            HlpSveImg(tFil,2);
        end
    end
    
    % every 100'th iteration, update trace plots, diagnostic plots and parameter plots
    if rem(iI,100)==0  
        
        % obtain trace plots
        figure(1); clf;
        PltTrc(C,iI,iI);
        tFil=sprintf('%s-TrcPrg%.2f',X.DatNam,100*X.Nep);
        HlpSveImg(tFil,2);
        
        % diagnostic plots
        figure(2); clf;
        PltDgn(C,iI);
        
        % plot parameters with distance and direction
        PltPrm(C,iI,NmbToPlt);
        
        % save end point for starting subsequent runs
        PrmPrg=C.Prm(iI,:)';
        NgtPrg=C.Ngt(iI,:)';
        ItrPrg=iI;
        tFil=sprintf('%s-PrmPrg%.2f.mat',X.DatNam,100*X.Nep);
        save(tFil,'PrmPrg','NgtPrg','ItrPrg');
        
        % print progress to terminal 
        fprintf(1,'Completed iteration %g. %g iterations in %g seconds\n',iI,100,toc);
        tic;
        
    end
    
    %on every 1000'th iteration, save parameter chain + associated trace
    %and diagnostic plots
    if rem(iI,1000)==0 || iI==mI+nI 
        
        % save whole MCMC chain
        DatNam=X.DatNam;
        Nep=X.Nep;
        tFil=sprintf('%s-Prm%.2f.mat',DatNam,100*Nep);
        save(tFil,'C','mI','nI','NmbToPlt','DatNam','Nep');
        
        % trace plots
        figure(1);
        clf; PltTrc(C,iI,iI); % save full trace plot
        HlpSveImg(sprintf('%s-trace%g%s',DatNam,100*Nep,datestr(now,30)),2);
        clf; PltTrc(C,iI,NmbToPlt); % save end of trace plot
        HlpSveImg(sprintf('%s-trace_end%g%s',DatNam,100*Nep,datestr(now,30)),2);
        
        % diagnostic plots
        figure(2);
        clf; PltDgn(C,iI);
        HlpSveImg(sprintf('%s-Dgn%g%s',DatNam,100*Nep,datestr(now,30)),2);    
        tStr=cell(3,1);
        tStr{1}=sprintf('%s-Prm%g%s',DatNam,100*Nep,datestr(now,30));
        tStr{2}=sprintf('%s-DrcAlp%g%s',DatNam,100*Nep,datestr(now,30));
        tStr{3}=sprintf('%s-DrcBet%g%s',DatNam,100*Nep,datestr(now,30));
        PltPrm(C,iI,NmbToPlt,tStr);
        
    end
    
end

return;