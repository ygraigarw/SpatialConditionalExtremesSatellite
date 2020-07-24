function HlpSveImg(Enw0,Fmt);
%function HlpSveImg(Enw0,Fmt);
%
%Philip Jonathan, Statistics & Chemometrics, Thornton. 
%"p*"-utilities 20091110
%
%Produce image from MATLAB figure window in TIFF and JPEG formats
%
%Input
% Enw0   character       Name of JPEG file
% Fmt  optional integer =1 => produce TIFF 
% Fmt  optional integer =2 => produce PNG

%
%Example
%HlpSveImg('Try') will generate "Try.jpg" from current MATLAB figure window
%Specify Fmt=1, and "Try.tif" will also be generated
%
%History
%20100912 - Refined header
if ~strcmp(get(gcf,'WindowStyle'),'docked')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end
set(gcf,'Color',[1 1 1]); %white background

% by default jpeg
Stn='jpeg';

if nargin==0;
	Enw0	= sprintf('GI_%s', datestr(now,30));    
end;

if nargin<=1;
    Fmt=0;
end
if nargin==2;
	if Fmt==1;
		Stn='tiff';
	elseif Fmt==2;
		Stn='png';
	elseif Fmt==0;
		Stn='jpeg';
    elseif Fmt==3;
        Stn='jpeg';
	end;
end

tStr=sprintf('Enw=''%s.%s'';',Enw0,Stn);eval(tStr);

try
    style = hgexport('factorystyle');
    hgexport(gcf, Enw, style, 'Format', Stn)
catch
    fprintf('Warning Picture not saved \n')
end

if Fmt==3;
saveas(gca,[Enw0,'.fig'],'fig')
end
%Normal completion

