function CQCBndL=DgnCQC(X);
%% function CQCBndL=DgnCQC(X);
%
% Conditional spatial extremes with delta-Laplace residuals:  ...
% Calculate conditional quantile constraint boundaries (Keef et al, 2013)
%
% Philip Jonathan, Rob Shooter, Emma Ross
%
%%

Nep=X.Nep;
DatNam=X.DatNam;
nLags=size(X.XRC,2);
%initialise empty cell array: each element will contain a boundary (alpha-beta coordinates) 
CQCBndL=cell(nLags,1); 


for iL=1:nLags  %loop over remote locations
    
    x=X.X0C;
    y=X.XRC(:,iL);
    
    n=200; 
    isOK=zeros(n,n);
    alp=linspace(-1,1,n);
    bet=linspace(-1,1,n);
    
    setGood=[];
    for iA=1:n
        for iB=1:n
            isOK(iA,iB)=CQC(alp(iA),bet(iB),[],[],x,y);
            if isOK(iA,iB)==1
                setGood=[setGood;[alp(iA) bet(iB)]];
            end
        end
    end
    
    isBnd=convhull(setGood(:,1),setGood(:,2));
    Bnd=[setGood(isBnd,1) setGood(isBnd,2)];
    Bnd=[1 0; Bnd(Bnd(:,1)>-0.5 & Bnd(:,2)>0,:)];
    xGrd=(0:0.01:1)';
    yGrd=interp1(Bnd(:,1),Bnd(:,2),xGrd);
    
    subplot(2,2,1);
    plot(x,y,'k.');
    xlabel 'Conditioning variate (Laplace scale)';
    ylabel 'Conditioned variate (Laplace scale)';
    title 'Sample';
    
    subplot(2,2,2);
    plot(xGrd,yGrd,'k.');
    xlabel '\alpha';
    ylabel '\beta';
    title 'CQC boundary';
    
    subplot(2,2,3);
    imagesc(alp,bet,isOK');
    xlabel('\alpha');
    ylabel('\beta');
    title 'Allowed domain (all pairs with given lag)';
    axis xy;
    
    subplot(2,2,4);
    imagesc(alp,bet,isOK');
    xlabel('\alpha');
    ylabel('\beta');
    axis xy;
    title 'Reduced domain';
    set(gca,'xlim',[0 1],'ylim',[0 1]);
    
    drawnow;
    
    %tStr=sprintf('%s-nep%g-lag%g',DatNam,Nep,iL);
    %pDatStm(tStr);
    HlpSveImg(sprintf('%s-CQCBnd-nep%g-lag%g',DatNam,round(Nep*100),iL),2);
    
    %save CQC boundary to output
    CQCBndL{iL}=[xGrd yGrd];
    
end
