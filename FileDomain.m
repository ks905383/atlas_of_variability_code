function VarDomain = FileDomain(filevar)
% FILEDOMAIN     get folder name for saving conventions.
%   VarDomain = FILEDOMAIN(filevar) gives [VarDomain] - the name of the
%   folder in which final figures and results are stored for saving
%   purposes - for the variable given by the string [filevar] (using CMIP5
%   variable abbreviations, found in [code_dir]/Varnames.csv). 
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%   
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 08/26/2016 (just added header, otherwise hasn't been
%   touched in geological ages)
 
if strcmp(filevar(1),'c')==1
    VarDomain = 'Clouds';
elseif strcmp(filevar(1),'r')==1
    VarDomain = 'Radiation';
elseif strcmp(filevar(1),'t')==1
    VarDomain = 'Temperature';
elseif strcmp(filevar(1),'m')==1
    VarDomain = 'SoilMoisture';
elseif strcmp(filevar(1),'u')==1 || strcmp(filevar(1),'v')==1
    VarDomain = 'Winds';
elseif strcmp(filevar(1),'h')==1
    VarDomain = 'HeatFlux';
elseif strcmp(filevar(1),'p')==1
    VarDomain = 'Precipitation';
elseif strcmp(filevar(1),'s')==1
    VarDomain = 'Sea';
end
end