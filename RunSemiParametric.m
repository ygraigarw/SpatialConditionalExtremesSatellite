%% CONTROL FILE Non Parametric -Conditional Spatial Extremes (NP CSE)
% 'Run-script' for Non-parametric conditional spatial extremes with delta-Laplace residuals
% Philip Jonathan, Rob Shooter, Emma Ross

% Script to generate sample data (or read data file) and make MCMC inference

%% STEPS to runs analysis
%1. Read the user guide which comes with this code
%2. Enter required input values/settings in the section bookended by "START USER INPUT / END USER INPUT"
%3. Move to 'Review Starting Solution' section to ensure a suitable
%starting solution has been chosen
%4. Run script (F5)

%% Specify new random seed for analysis
rand('state',sum(100*clock)); %new starting point for random numbers based on clock

%% Set path 
PthThisFil = matlab.desktop.editor.getActiveFilename;  %full pathname of file currently open in users console (wherever this function GetPthActvFil is being run from)  
PosDlm = strfind(PthThisFil,'\');  %find delimiters (/) in directory
HomPth= extractBefore(PthThisFil,PosDlm(length(PosDlm))); % return directory in which the active script sits

%% Input arguments

%START USER INPUT >>

%Input data file
X.DatFil=fullfile(HomPth,'Satellite_Preprocessed_NS_SWNE'); %input data file to use (whitened Hs data?)

% X: peaks data-related
X.Nep=0.8; %non-exceedance probability for CSE model
X.CndLct=1; %index of conditioning location
% C: spatial conditional extremes model related
C.IsCQC=1; %=1 == Impose conditional quantile constraints; =0 == No constraints
C.StrAgn=0; %=1 == Start the analysis from scratch; =0 == Do not restart the analysis, i.e. use previously-saved progress
C.nI=10; %number of MCMC iterations
C.Rpl=1; %replicate number (to distinguish replicate analyses)
C.nH=10; %number of distance nodes
C.nA=1; %number of directional nodes
C.HMxm=12; %maximum interlocation distance; in units of EarthRadius*pi/180; one unit is 111.2km

Scl1=100; 
Scl2=2;
% Starting solution
C.Prm0=[...
    ones(C.nH*C.nA,1)*0.5;...  %for all alphas and angles
    ones(C.nH*C.nA,1)*0.1;... %for all betas and angles
    ones(C.nH,1)*0.3; %for all mus
    ones(C.nH,1)*1.0; %for all sigma
    ones(C.nH,1)*1.5; %for all deltas
    0.6;0.56;0.5];  %rho1, rho2, (rho3 - could take out)

StartingSolutionOk = true; %set to true when have found a starting solution C.Prm0 which satisfies CQC boundaries ...
%...(if CQC boundaries turned off, ie. C.IsCQC=0; doesn't matter what this
%is set to)

% OPTIONAL SETTINGS-ADJUSTMENTS 
% Parameter limits
C.Lmt.A=[-0.1 1]; %limit for each distance and direction for alpha SCE parameter
C.Lmt.B=[-0.1 1]; %limit for each distance and direction for beta SCE parameter
C.Lmt.M=[-1 1]; %limit for each distance for mu SCE parameter
C.Lmt.S=[0 sqrt(2)+0.1]; %limit for each distance for sigma SCE parameter
C.Lmt.D=[0.5 2.5]; %limit for each distance for delta SCE parameter
C.Lmt.R=[0 1]; %limit for the parameter rho in the residual process
% MCMC Setings
C.NgtStr=0.1; %candidate random walk standard deviation
C.NmbToPlt=1000; %number of iterations (at end of the chain) to use for plotting
C.AdpItr1=500; %(adaptive MCMC setting) no. iterations for fixed nugget. Warning: must take minimum value of 100 
C.AdpItr2=500; %(adaptive MCMC setting) no. iterations for adaptive nugget. Warning: must take minimum value of 100 
C.AdpBet=0.05; %(adaptive MCMC setting) proportion of proposal distribution made up of adaptive nugget 
%parameters used in GetRsdCrr function, suitable for satellite application
C.Scl1=100;    
C.Scl2=2;    

% << END USER INPUT

C.Dlt=C.HMxm/(C.nH-1); %distance increment

%% Folder/Directory Handling

%add code to path
CodPth=fullfile(HomPth,'Code');
addpath(CodPth);

%Create directory for analysis and move to this folder 
X.DatNam=sprintf('Nep%gnH%gnA%gCl%gCqc%gRpt%g',100*X.Nep,C.nH,C.nA,X.CndLct,C.IsCQC,C.Rpl); %define name for analysis
X.AnlPth=fullfile(HomPth,X.DatNam); %define analysis path
if exist(X.AnlPth,'dir')==0 %create folder if it doesnt exist
    mkdir(X.AnlPth);
end
cd(X.AnlPth); %change directory, move to analysis folder

fprintf(1,'Spatial conditional extremes %s started at %s\n',X.DatNam,datestr(now,30));  %print to matlab terminal: time job started

%% Format Peaks Data 
% A sample datafile containing structure 'Pk' is provided. To load your own
% data into the required format: 
% * Populate Dat.Unf [n x p array] with n observations at p locations (pre-'whitened' peak Hs data, on uniform scale)
% * Popualate Dat.Lct [p x 2 array] with longitude and latitude for each of the p locations

%START OPTIONAL USER INPUT: INPUT PEAKS DATA >>
load(X.DatFil,'Pk'); %load data file
Dat.Unf=Pk.HsUnf; %peaks-data on uniform marginal scale 
Dat.Lct=Pk.Lct; %Longitudes and latitudes of locations for analysis
clear Pk;
% << END OPTIONAL USER INPUT

% Make sure conditioning location is first location
AllLct=(1:1:size(Dat.Unf,2))'; %specify all locations to consider (e.g. if want to use only every other location, AllLct=(1:2:size(Dat.Unf,2))
IsCnd=ismember(AllLct,X.CndLct); %confirm that conditioning location occurs in the data file
X.LctNmb=[X.CndLct;AllLct(IsCnd==0)]; %rearrange locations, conditioning location first
DatU=Dat.Unf(:,X.LctNmb); %data should be on uniform marginal scale (obtained by fitting a marginal model to peaks data and applying the probability integral transform). 
X.L0=Dat.Lct(X.LctNmb,:); %p locations on the x-axis in 2D space

% Set non-exceedance threshold for conditioning location
X.ThrLpl=-sign(X.Nep-0.5).*log(1 - 2.*abs(X.Nep-0.5)); %non-exceedance threshold on Laplace-scale
X.DatL=-sign(DatU-0.5).*log(1 - 2.*abs(DatU-0.5)); %Uniform to Laplace
X.p=size(X.DatL,2); %total number of locations
X.q=X.p-1; %number of remote locations

% Co-ordinates of reference and remote locations
X.Rfr=X.L0(1,:); %reference locations
X.Rmt=X.L0((2:X.p),:); %remote locations

% Data (laplace scale) for reference and remote locations
tKep=X.DatL(:,1)>X.ThrLpl; %indicator for threshold exceedances at reference location
X.X0C=X.DatL(tKep,1); %threshold exceedances at reference location
X.XRC=X.DatL(tKep,(2:X.p)); %threshold exceedances at remote location


%% Plot data locations
% Reference image, where on the transect does data indexed by XX appear
set(groot,'DefaultAxesColorOrder',jet(X.q));
figure(1); clf; hold on;
plot(X.L0(:,1),X.L0(:,2),'w.');
xlabel('Longitude')
ylabel('Latitude')
for j=1:X.p
    text(X.L0(j,1),X.L0(j,2),sprintf('%g',X.LctNmb(j)));
end
HlpSveImg(sprintf('%s-Locations%g%s',X.DatNam,100*X.Nep,datestr(now,30)),2);

%% Review starting solution
% Need to be careful with starting solution when CQCs imposed, since starting solution can violate CQCs
% Plotting CQC boundaries helps identify suitable values


if C.IsCQC==1 % if have conditional quantile constraints switched on
    
    % --- Load or define CQC boundaries
    tFil=sprintf('%s-CQCBnd.mat',X.DatNam); %filename for conditional quantile constraints
    if exist(fullfile(cd,tFil),'file')==2  %if file already exists and don't want to regenerate it then load
        fprintf(1,'Using existing CQC boundaries\n');
        load(tFil,'CQCBndL','CQCNepLst');    
    else % if don't have CQC boundaries already or want to regenerate, calculate ...
        %...and check that your starting value satifies them
        fprintf(1,'Generating new CQC boundaries\n');
        CQCBndL=DgnCQC(X); %plot the domains of CQC constraints for lagDatL
        CQCNepLst=X.Nep;
        save(tFil,'CQCBndL','CQCNepLst');
        fprintf(1,'Saved CQC boundaries\n');
    end
    C.CQCBndL=CQCBndL;
    
    % ---- Check starting value: once have landed on C.Prm0 which is
    % feasible, set above StartingSolutionOk=true and code will move onto
    % MCMC
    if ~StartingSolutionOk 
        fprintf(1,'Action: Check starting solution using figure, when starting solution feasible set StartingSolutionOk to TRUE and rerun to continue analysis');
        tClr = lines(X.q); 
        figure;hold on;
        for iL=1:X.q %loop over remote locations, plotting keef boundaries
            if C.IsCQC==1
                plot(C.CQCBndL{iL}(:,1),C.CQCBndL{iL}(:,2),'color',tClr(iL,:));
            end
        end
        stAlps = C.Prm0(1:C.nH*C.nA); %starting alphas
        stBets = C.Prm0(C.nH*C.nA+1:2*C.nH*C.nA);% starting betas
        plot(stAlps,stBets,'.k', 'MarkerSize',10); %plot starting parameters
        HlpFnt; HlpAxsLmt;
        xlabel('\alpha')
        ylabel('\beta')
        title('Check starting value satisfies CQC: should be below-left of curves')
        
        %stop code here whilst finding suitable starting value C.Prm0. To
        %move onto MCMC fitting, set StartingSolutionOk=true. 
        return  

    end
end      



%% Run MCMC code
MCMC(X,C);

%% Return to home folder
cd(HomPth);

return;