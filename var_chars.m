function varargout = var_chars(var_idx,varargin)
% VAR_CHARS     get information about a variable
%   [filevar,filevarFN,freq] = VAR_CHARS(var_idx) gives the CMIP5
%   designated / original variable name [filevar], the internally
%   designated (post-processed) variable name [filevarFN], and the
%   frequency [freq] of the climate variable described by [var_idx] (if
%   [var_idx] is an integer, this corresponds to the [var_idx]th row in
%   [code_dir]/Varnames.csv. If [var_idx] is a string, VAR_CHARS
%   will find the row in Varnames.csv whose second column (post-processed
%   variable name) matches [var_idx]. If more than one frequency is found
%   for the string-set variable, a warning is given, and the first (by row
%   order in Varnames.csv) is used).
%
%   [filevar,filevarFN,freq,vardesc,units,varunits,data_type,freqdesc,frequnit]
%   = VAR_CHARS(var_idx) additionally outputs a short-form variable
%   decsription [vardesc], the units [units] it is stored in, the units
%   squared [varunits], the data_type [data_type], a short-form frequency
%   description, and the number of data points [frequnit] in a year.
%
%   [__] = VAR_CHARS(filevarFN,freq) allows finding variable
%   characteristics with a specified (string) frequency as well. 
%
%   All outputs except [frequnit] are strings. 
%   [data_type] gives the following:
%       'percent'   - percentage data
%       'temp'      - temperature data
%       'ratio'     - ratio data not covered by 'percent' or 'temp'
%       'interval'  - interval data
%   These [data_type] outputs are used in determining colormaps and
%   eligibility for coefficient of variation processing (see also
%   Maps_pi_extras and Maps_pi_flat)
%
%   Example:
%   [___] = VAR_CHARS(var_idx) with var_idx = 11 or 'tas' or 'tas_day' or
%   'tas','day' will return
%       filevar = 'tas'
%       filevarFN = 'tas'
%       freq = '_day_'
%       vardesc = 'Near-Surface Air Temperature'
%       units = '( K )'
%       varunits = '( K^2 )'
%       data_type = 'temp'
%       freqdesc = 'Daily '
%       frequnit = 365
%
%   Comments: the difference between filevar and filevarFN lies in
%   different subsetting of variables between CMIP5 and the internal file
%   system. For example, 850mb meridional wind is stored as one of several 
%   pressure elevation bands in the CMIP5 format in the variable (filevar)
%   "va", while it is pre-processed in this filesystem into an independent
%   variable "va850". For more info, see the program Saves.m (link below).
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also NAME_CHARS, SAVES
%   
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified: 10/14/2016

%Set regex delimiters
delim1 = '\_';

%Load file with variable characteristics
various_defaults = matfile('various_defaults.mat');
VarnamesFile = fopen([various_defaults.code_dir,'Varnames.csv']);
if VarnamesFile<0
    error('VAR_CHARS:NoVarnames',['Could not find the file Varnames.csv in ',...
        various_defaults.code_dir,'. Please ensure this file exists at this location.'])
else
    Varnames = textscan(VarnamesFile,'%s','Delimiter',',');
end
fclose(VarnamesFile);
Varnames = reshape(Varnames{1},[8 (length(Varnames{1})/8)]);

%Allow for string identifiers of var_idx (not as robust)
if isa(var_idx,'char')
    %Support for input of the form '[var]_[freq](_)'
    if ~isempty(strfind(var_idx,'_'));
        idx_tmp = regexp(var_idx,delim1,'split');
        var_idx = idx_tmp{1};
        if isempty(varargin)
            varargin{1} = idx_tmp{2};
        end
    end
    %Find desired variable in relevant Varnames column (2, for filevarFN,
    %aka post-processed filename (for more info, see Saves))
    var_idx_tmp = find(cellfun(@(x) strcmp(var_idx,x),Varnames(2,:)));
    %If string found in Varnames, continue, else throw error
    if ~isempty(var_idx_tmp)
        %Can insert frequency as a second input to refine between two
        %possible variable ids
        if ~isempty(varargin)
            freq_idx = find(~cellfun(@isempty,strfind(Varnames(3,:),varargin{1})));
            var_idx = intersect(freq_idx,var_idx_tmp); var_idx = var_idx(1);
        else %if no frequency specified, just pick first one
            if length(var_idx_tmp)>1
                warning('var_chars:MultOpts',['There are multiple possible data frequencies listed for "',...
                    var_idx,'," the first one in the order of Varnames.csv (row ',num2str(var_idx_tmp(1)),') will be used'])
                var_idx = var_idx_tmp(1);
            else
                var_idx = var_idx_tmp;
            end
        end
    else
        error('var_chars:InvalidInput',['"',var_idx,'" is not a supported variable index string'])
    end
end

%Define file variable identifiers
filevar=Varnames{1,var_idx};
filevarFN=Varnames{2,var_idx};
freqStr = regexp(Varnames{3,var_idx},delim1,'split');
freq=strcat('_',freqStr{2},'_'); clear freqStr;

vardesc = Varnames{4,var_idx};
units = Varnames{5,var_idx};
varunits = Varnames{6,var_idx};
data_type = Varnames{7,var_idx};

if regexp(freq,'.[Dd]ay.')
	freqdesc='Daily ';
	frequnit=365;
elseif regexp(freq,'.[Mm]on.')
	freqdesc='Monthly ';
	frequnit=12;
elseif regexp(freq,'.yrr.')
    freqdesc = 'Yearly ';
    frequnit = 1;
elseif regexp(freq,'.3hr.')
    freqdesc = '3-hour ';
    frequnit = 2920;
elseif regexp(freq,'.6hr.')
    freqdesc = '6-hour ';
    frequnit = 1460;
elseif regexp(freq,'.12hr.')
    freqdesc = '12-hour ';
    frequnit = 730;
end

varargout{1} = filevar; varargout{2} = filevarFN; varargout{3} = freq;
if nargout > 3
    varargout{4} = vardesc; varargout{5} = units; varargout{6} = varunits;
    varargout{7} = data_type; varargout{8} = freqdesc; varargout{9} = frequnit;
end

end