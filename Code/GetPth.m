function HomPth=GetPthActvFil()
%Retrieve directory location of the matlab file active in matlab console (tab the user is viewing)

PthThisFil = matlab.desktop.editor.getActiveFilename;  %full pathname of file currently open in users console (wherever this function GetPthActvFil is being run from)  
PosDlm = strfind(PthThisFil,'\');  %find delimiters (/) in directory
HomPth= extractBefore(PthThisFil,PosDlm(length(PosDlm))); % return directory in which the active script sits

end