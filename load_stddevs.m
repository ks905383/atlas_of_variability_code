function varargout = load_stddevs(var_idx,model,exp,varargin)
% LOAD_STDDEVS    load variability values from files without worrying about 
%                 source file convention.
%
%   StdDevs = LOAD_STDDEVS(var_idx,model,exp) outputs the variability
%   values of the the variable defined by variable index or name (see
%   var_chars documentation for more details on possible var_idx inputs)
%   [var_idx], model [model], and experiment [exp]. [StdDevs] is a (number
%   of frequency bands) x 1 struct with fields .Name (for graph titles),
%   .ID (shortform), .Size (period start and end points for each band; [0
%   Inf] or [0 Inf] for all frequencies for monthly, daily data
%   respectively), and .Data (containing the data). These are identical to
%   the outputs of Variability and Variability_seasons and the new (as of
%   date last modified below) filesystem and structure (and are simply
%   outputted if the origin file follows the new file structure, or created
%   if not). Supported file conventions include StdDevs, [ID]avgDev, and
%   [ID]avgDev[exp] (all except StdDevs being depreciated in newer code).
%
%   StdDevs = LOAD_STDDEVS(var_idx,model,exp,season) instead outputs the
%   same, but limited to the season defined by [season] (by default,
%   [season] = 'all'. Season IDs have to be all uppercase, but are
%   otherwise just set by the relevant position of the file extension, i.e.
%   'DJF','JJA','MAMJJA', etc.).
%
%   [StdDevs,lat,lon] = LOAD_STDDEVS(var_idx,model,exp,...) also outputs
%   the latitude and longitude vectors of the model and variable.
%
%   [lat,lon] = LOAD_STDDEVS(var_idx,model,exp,...) only outputs the
%   latitude and longitude vectors of the model and variable (slightly
%   faster when the full variability struct isn't needed). 
%
%   property = LOAD_STDDEVS(___,'property',[prop_name]) instead outputs
%   a cell containing the values of the struct field given by
%   [prop_name] - i.e. LOAD_STDDEVS(___,'property','ID') will return, i.e., 
%
%       property = {'HF','MF','LF','XF','seas_cycle','full_noDS','mad',...}
%
%   if those are the frequency band ids contained in the file to be read.
%
%   filename = LOAD_STDDEVS(___,'filename') instead outputs a string giving
%   the filename that would be loaded.
%
%   Sample use:
%   [StdDevs,lat,lon] = LOAD_STDDEVS(11,'NorESM-1M','rcp85','DJF') outputs
%   a struct [StdDevs] containing each calculated variability parameter
%   (usually standard deviation split up by frequency band contributions
%   plus a few auxiliary calculations such as the strength of the seasonal
%   cycle, etc.) at each location of daily reference height temperature 
%   (var_idx = 11, or alternatively 'tas' or 'tas_day') in the RCP8.5
%   scenario for winter, over the entirety of the time series saved in this 
%   file system (usually 30 years/season), in addition to two vectors [lat]
%   and [lon] containing the lat and lon values of each column and row in
%   the grid.
%
%   Notes: - by default, the first use_convention will be used.
%          - warnings are thrown if more than one variable matching
%            '*StdDev*' are found, with those variables listed (shouldn't
%            happen though, tbh).
%          - by default, name_chars warnings about multiple files for a 
%            given variable-model-experiment combination are turned off.
%          - for general coding flexibility, [season_folder] (option
%            below) does not change running of function if [season] = 'all'
%
%   LOAD_STDDEVS(...,'[flag]',[params],...) modify program run as below:
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
%   See also LOAD_MEANS, NAME_CHARS,
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 11/26/2017

%Set defaults
various_defaults = matfile('various_defaults.mat');
folder_season = various_defaults.season_dir;
use_convention = 1;
strtyr_find = []; endyr_find = [];
keep_warnings = false;
struct_pop_values = {'HF','MF','LF','XF','Full';...
    5,30,365,36500,0;...
    NaN,NaN,12,1200,0;...
    '< 5 days','5 - 30 days','30 - 365 days','> 365 days','all frequencies';...
    NaN,NaN,'< 1 year','> 1 year','all frequencies'};
property_out = [];
output_filename = false;

%Set season 'all' to be the default version. Based on no other input flags
%being 'all' or all upper case, but I think that's pretty easy to keep. 
if ~isempty(varargin) && isa(varargin{1},'char') && (strcmp(varargin{1},'all') || all(isstrprop(varargin{1},'upper')))
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
            case {'end_year'} %NOT YET DEVELOPED
                endyr_find = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'run'} %NOT YET DEVELOPED
                run_set = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'keep_warnings'}
                keep_warnings = true;
            case {'property'}
                property_out = varargin{in_idx+1};
            case {'filename'}
                output_filename = true;
        end
    end
end

%Set warnings from name_chars to on or off
if keep_warnings; warning('on','name_chars:MultConvs');
else warning('off','name_chars:MultConvs'); %Unnecessary if you're clear on what you're doing
end
warning('off','name_chars:FileNotFound'); %This is covered by the error below

%Get basic file information
[~,filevar,freq] = var_chars(var_idx);
%Get file locations and characteristics
if strcmp(season,'all') %Get location standards for non-seasonal data
    folder = '/';
    filename_add_open = [];
    [~,run,strtyr,endyr,num_conventions] = name_chars(model,{exp},filevar,freq,'use_convention',use_convention);
    if ~isempty(strtyr_find) %If explicitly setting start year
        if strtyr{1} ~= strtyr_find
            if num_conventions > 1 %Find right file convention
                n = use_convention+1;
                while strtyr{1} ~= strtyr_find
                    if n > num_conventions; error('LOAD_STDDEVS:FileNotFound',['No file with a start year of ',num2str(strtyr_find),' was found. No variability data was loaded.']); end
                    [~,run,strtyr,endyr,num_conventions] = name_chars(model,{exp},filevar,freq,'use_convention',n);
                    n = n+1;
                end
            else
                error('LOAD_STDDEVS:FileNotFound',['No file with a start year of ',num2str(strtyr_find),' was found. No variability data was loaded.']);
            end
        end
    end
else %Get location standards for seasonal data
    filename_add_open = ['_',season];
    folder = folder_season; dir_seas = [various_defaults.proc_data_dir,model,folder];
    [~,run,strtyr,endyr,num_conventions] = name_chars(model,{exp},filevar,freq,'directory',dir_seas,'file_end',['*Power_',season,'.mat'],'use_convention',use_convention);
    if ~isempty(strtyr_find) %If explicitly setting start year
        if strtyr{1} ~= strtyr_find
            if num_conventions > 1 %Find right file convention
                n = use_convention+1;
                while strtyr{1} ~= strtyr_find
                    if n > num_conventions; error('LOAD_STDDEVS:FileNotFound',['No file with a start year of ',num2str(strtyr_find),' was found. No variability data was loaded.']); end
                    [~,run,strtyr,endyr,num_conventions] = name_chars(model,{exp},filevar,freq,'directory',dir_seas,'file_end',['*Power_',season,'.mat'],'use_convention',n);
                    n = n+1;
                end
            else
                error('LOAD_STDDEVS:FileNotFound',['No file with a start year of ',num2str(strtyr_find),' was found for ',filevar,freq,model,'_',exp,'. No variability data was loaded.']);
            end
        end
    end
end

%Load file
filename = strcat(various_defaults.proc_data_dir,model,folder,filevar,freq,model,'_',exp,'_',run{1},num2str(strtyr{1}),'-',num2str(endyr{1}),'_sqrtPower',filename_add_open,'.mat');
if exist(filename,'file')
    file1 = matfile(filename);
else
    error('LOAD_STDDEVS:FileNotFound',['The file ',filename,' was not found. No variability data was loaded.'])
end

if (nargout ~= 2 && ~output_filename) || ~isempty(property_out)
    load(filename)
    %Deal with different file conventions to extract data
    if ~exist('StdDevs','var')
        StdDevs_tmp = who('*StdDev*');
        if isempty(StdDevs_tmp)
            StdDevs_tmp = who('*avgDev*');
            if ~isempty(StdDevs_tmp)
                StdDevs(length(StdDevs_tmp)).Name = 'Initialize'; %Initialize struct in memory
                for var = 1:length(StdDevs_tmp)
                    StdDevs(var).ID = strtok(StdDevs_tmp{var},'avg');
                    StdDevs(var).Data = eval(StdDevs_tmp{var});
                    pop_idx =  find(strcmp(StdDevs(var).ID,struct_pop_values(1,:)));
                    if strcmp(freq,'_day_')
                        StdDevs(var).Size = struct_pop_values{2,pop_idx};
                        StdDevs(var).Name = struct_pop_values{4,pop_idx};
                    else
                        StdDevs(var).Size = struct_pop_values{3,pop_idx};
                        StdDevs(var).Name = struct_pop_values{5,pop_idx};
                    end
                end
            else
                error('LOAD_STDDEVS:StdNotFound',['No variable with the name "*StdDev*" or "*avgDev*" found in ',filename,' please check ',filename,' for filesystem compatibility and existence of data'])
            end
        elseif length(StdDevs_tmp)>1
            warning('LOAD_STDDEVS:TooManyStds',['More than one variable with the name "*StdDev*" found, please check ',filename,' for filesystem compatibility. The variables found are: ',sprintf('\n%s',StdDevs_tmp{:})])
        end
    end
    
    if ~isempty(property_out)
        
        prop_eval = ['{StdDevs(:).',property_out,'}'];
        varargout{1} = eval(prop_eval);
    else
        varargout{1} = StdDevs;
        
        if nargout > 1
            varargout{2} = lat;
            varargout{3} = lon;
        end
    end
else
    
    if output_filename
        varargout{1} = filename;
    else
        varargout{1} = file1.lat;
        varargout{2} = file1.lon;
    end
end


%Reset warning messages
warning('on','name_chars:MultConvs'); 
warning('on','name_chars:FileNotFound');
end




