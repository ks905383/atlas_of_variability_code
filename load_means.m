function varargout = load_means(var_idx,model,exp,varargin)
% LOAD_MEANS    load mean values from files without worrying about source
%               file convention.
%
%   LocalMean = LOAD_MEANS(var_idx,model,exp) outputs the local mean
%   value of the variable defined by variable index or name (see var_chars
%   documentation for more details on possible var_idx inputs) [var_idx],
%   model [model], experiment [exp]. This function supports LocalMeans,
%   Means[exp], LocalMean, and any other file convention that has saved the
%   [LocalMean] data as a variable containing the string '*Mean*'.
%
%   LocalMean = LOAD_MEANS(var_idx,model,exp,season) instead outputs the
%   local mean value of the variable in the season defined by [season] (by
%   default, [season] = 'all'. Season IDs have to be all uppercase, but are
%   otherwise just set by the relevant position of the file extension, i.e.
%   'DJF','JJA','MAMJJA', etc.). 
%
%   [LocalMean,lat,lon] = LOAD_MEANS(var_idx,model,exp,...) also outputs
%   the latitude and longitude vectors of the model and variable.
%
%   Notes: - by default, the first use_convention will be used.
%          - warnings are thrown if more than one variable matching
%            '*Mean*' are found, with those variables listed (usually
%            happens in files that had the data saved twice to be backwards
%            compatiable with the old MeansPI/etc. filesystem, in which
%            case the warning can safely be ignored).
%          - by default, name_chars warnings about multiple files for a 
%            given variable-model-experiment combination are turned off.
%          - for general coding flexibility, [season_folder] (option
%            below) does not change running of function of [season] = 'all'
%
%   Sample use:
%   means = LOAD_MEANS(11,'NorESM1-M','rcp85','all'); gives the local mean
%   value of daily reference height temperature (var_idx = 11, or 
%   alternatively 'tas' or 'tas_day') in the RCP8.5 scenario over the
%   entirety of the time series saved in this file system (usually 30
%   years)
%
%   LOAD_MEANS(...,'[flag]',[params],...) modify program run as below:
%       'season_folder',[char]      - changes folder in
%                                     [proc_data_dir]/[model] in which
%                                     seasonal data is looked for (def:
%                                     [seasonal_dir])
%       'use_convention',[#]        - specifies which timeframe of file to
%                                     look for by its spot in the
%                                     system-returned order in the target
%                                     directory (alphabetical in linux and,
%                                     honestly, hopefully any OS you're
%                                     using). Error thrown if greater than
%                                     number of conventions. Equivalent to
%                                     (and gets piped right into)
%                                     name_chars 'use_convention' flag.
%       'start_year',[#]            - specifies start year to look for.
%                                     Error thrown if not found. 
%       'keep_warnings'             - does not suppress the name_chars
%                                     multiple convention warnings (useful
%                                     if unclear if right timeframe is
%                                     being outputted)
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also LOAD_STDDEVS, NAME_CHARS,
%   
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 11/26/2017

%Set defaults
various_defaults = matfile('various_defaults.mat');
folder_season = various_defaults.season_dir;
use_convention = 1;
strtyr_find = [];
keep_warnings = false;

%Set season 'all' to be the default version. Based on no other input flags
%being 'all' or all upper case, but I think that's pretty easy to keep. 
if ~isempty(varargin) && isa(varargin{1},'character') && (strcmp(varargin{1},'all') || all(isstrprop(varargin{1},'upper')))
   season = varargin{1};
   varargin = varargin(2:end);
else
    season = 'all';
end

%Set optional flag behavior
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'season_folder'}
                folder_season = varargin{in_idx+1};
            case {'use_convention'}
                use_convention = varargin{in_idx+1};
            case {'start_year'}
                strtyr_find = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'keep_warnings'}
                keep_warnings = true;
        end
    end
end

%Set warnings from name_chars to on or off
if keep_warnings; warning('on','name_chars:MultConvs');
else warning('off','name_chars:MultConvs'); %Unnecessary if you're clear on what you're doing
end
warning('off','name_chars:FileNotFound'); %This is covered by the error below

%Get basic file information
[~,filevar,freq,~,~,~,~,~,~] = var_chars(var_idx);
%Get file locations and characteristics
if strcmp(season,'all') %Get location standards for non-seasonal data
    folder = '/';
    filename_add_open = [];
    mean_fileid = '_LocalMeans';
    [~,run,strtyr,endyr,num_conventions] = name_chars(model,{exp},filevar,freq,'use_convention',use_convention);
    if ~isempty(strtyr_find) %If explicitly setting start year
        if strtyr{1} ~= strtyr_find
            if num_conventions > 1 %Find right file convention
                n = use_convention+1;
                while strtyr{1} ~= strtyr_find
                    if n > num_conventions; error('load_means:FileNotFound',['No file with a start year of ',num2str(strtyr_find),' was found. No means were loaded.']); end
                    [~,run,strtyr,endyr,num_conventions] = name_chars(model,{exp},filevar,freq,'use_convention',n);
                    n = n+1;
                end
            else
                error('load_means:FileNotFound',['No file with a start year of ',num2str(strtyr_find),' was found. No means were loaded.']);
            end
        end
    end
else %Get location standards for seasonal data
    filename_add_open = ['_',season];
    %mean_fileid = '_InterannualStdDev';
    mean_fileid = '_LocalMeans';
    folder = folder_season; dir_seas = [various_defaults.proc_data_dir,model,folder];
    [~,run,strtyr,endyr,num_conventions] = name_chars(model,{exp},filevar,freq,'directory',dir_seas,'file_end',['*Power_',season,'.mat'],'use_convention',use_convention);
    if ~isempty(strtyr_find) %If explicitly setting start year
        if strtyr{1} ~= strtyr_find
            if num_conventions > 1 %Find right file convention
                n = use_convention+1;
                while strtyr{1} ~= strtyr_find
                    if n > num_conventions; error('load_means:FileNotFound',['No file with a start year of ',num2str(strtyr_find),' was found. No means were loaded.']); end
                    [~,run,strtyr,endyr,num_conventions] = name_chars(model,{exp},filevar,freq,'directory',dir_seas,'file_end',['*Power_',season,'.mat'],'use_convention',n);
                    n = n+1;
                end
            else
                error('load_means:FileNotFound',['No file with a start year of ',num2str(strtyr_find),' was found. No means were loaded.']);
            end
        end
    end
end

%Load file
filename = strcat(various_defaults.proc_data_dir,model,folder,filevar,freq,model,'_',exp,'_',run{1},num2str(strtyr{1}),'-',num2str(endyr{1}),mean_fileid,filename_add_open,'.mat');
if exist(filename,'file')
    load(filename); %Means
else
    error('load_means:FileNotFound',['The file ',filename,' was not found. No means were loaded.'])
end
%Deal with different file conventions to extract data
Means_tmp = who('*Mean*');
if isempty(Means_tmp)
    error('load_means:MeansNotFound',['No variable with the name "*Mean*" found in ',filename,' please check ',filename,' for filesystem compatibility and existence of data'])
elseif length(Means_tmp)>1
    warning('load_means:TooManyMeans',['More than one variable with the name "*Mean*" found, please check ',filename,' for filesystem compatibility. The variables found are: ',sprintf('\n%s',Means_tmp{:})])
end


varargout{1} = eval([Means_tmp{1},';']);

if nargout > 1
    varargout{2} = lat;
    varargout{3} = lon;
end

%Reset warning messages
warning('on','name_chars:MultConvs'); %Unnecessary if you're clear on what you're doing
warning('on','name_chars:FileNotFound'); %This is covered by the error below

end