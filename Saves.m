function Saves(VarIndices,modelArray,varargin)
% SAVES     Extract data from netCDF files (CMIP5 standard conventions)
%           into .mat files
%
%   SAVES(VarIndices,modelArray) extracts 31 years of data from NetCDF
%   files saved in [raw_data_dir]/[model]/ containing data for
%   each variable and model given by [VarIndices] and [modelArray]. 30
%   years are saved in [proc_data_dir]/[model]/, 31 years are
%   saved in [proc_data_dir]/[model]/[season_dir]/. The general
%   form of each saved file (for a 3-dimensional variable reported over
%   lat, lon, and time) includes a variable [Raw] with the nlon x nlat x
%   time time series of the desired variable, variables [lon] and [lat]
%   giving the longitude and latitude values for the variable, and
%   variables [lon_bnds] and [lat_bnds] giving the extent of the pixels in
%   for the geographic grid. 
%
%   (All directories mentioned above are set by [various_defaults.m])
%
%   SAVES searches for the primary variable given by the row index of
%   Varnames.csv inputted as [VarIndices]. If the primary variable name (or
%   'lat', 'lon', 'lat_bnds', or 'lon_bnds') is not found in the original
%   NetCDF file, the program run is interrupted and the available variables
%   are shown, allowing the user to manually pick which variable to assign
%   to 'Raw', 'lat', 'lon', 'lat_bnds', or 'lon_bnds', or to interrupt
%   (skip) the program run. This is to take into account file conventions
%   in which the variables were saved with incorrect or out of convention
%   filenames. If the difference in filename is known from the start, the
%   program can be set to specifically look for that different filename
%   (see flags below). If the file is known to not contain [lat_bnds] or
%   [lon_bnds], or they are not desired in the final output file, use the
%   flag 'skip_bounds' set to false (see flags below). 
%
%   LEAP YEARS: SAVES removes the 366th day from any leap year (by
%   gregorian or julian leap year standards, for corresponding leap year
%   calendars only)
%
%   SUPPORTED CALENDARS: 365_day, 360_day, gregorian, proleptic_gregorian,
%   no_leap, and julian (all CMIP5-used calendars, as detailed <a href =
%   "https://esgf.github.io/esgf-swt/data/2015/08/13/Do-all-CMIP5-models-use-the-same-calendar.html">here</a>)
%
%   SUPPORTED DATA FREQUENCIES: daily (including cfDay, etc.), monthly
%   (including Amon, Lmon, etc.), and yearly (supporting *ann* and *yrr*)
%
%   SUPPORTED/EXPECTED NETCDF FILENAME CONVENTION:
%       [filevar]_[freq]_[model]_[exp]_[run]_[startdate]_[enddate].nc
%   with [startdate] and [enddate] being Y{1-4} (yearly data), YYYYMM
%   (monthly data), or YYYYMMDD (daily data)
%
%   SAVES is designed to minimize memory usage in MATLAB, and by default
%   extracts the minimum required for the desired time series length from
%   the netCDF file(s). Total memory used will be slightly more than that
%   needed for an array of size ~ [(nlon+1) x (nlat+1) x (number of time
%   points to save)] (the +1s give an estimate of the size of auxiliary
%   arrays and other variables used in the saving process). Additionally,
%   variables are saved in order of size (lat, lon, Raw) to guard against
%   excessive overhead when loading single variables caused by
%   decompressing behavior in certain file formats (i.e. v6 and v7, in
%   which files are decompressed in order they are saved, even if only the
%   last is desired).
%
%   Function requirements (on path): var_chars
%
%   SAVES(...,'[flag]',[params],...) modify program run as below:
%       'replace_files',[log]   - choose whether to replace existing files
%                                 with the same filename as outputted by
%                                 SAVES
%       'experiments',[cell]    - manually set experiments to process
%                                 (default {'rcp85','historical',
%                                 'piControl'})
%       'num_years',[int],([int])- how many years of data to take from the
%                                 end of the group of netcdf files and save
%                                 in a _DATA.mat file, by default 30.
%                                 Second integer is how many years to save
%                                 into the seasonal directory; if left
%                                 empty, the additional file is not
%                                 created.
%       'year_range',[int],[int]- manually choose start and end years for
%                                 the file to be saved (i.e. 1999,2005) (no
%                                 additional file is saved in the seasonal
%                                 directory in this case/save_case is set
%                                 to 'std')
%       'save_case',[char]      - choose to only save standard saves or
%                                 only alt (seasons by def) saves, through
%                                 'std', 'alt', or 'all'
%       'skip_bounds',[log]     - whether to save/look for lat/lon bounds
%                                 in the netcdf file (def: true) - use if
%                                 the netcdf files don't have pixel bounds
%                                 saved
%       'testing'               - adds '_TEST' to filename (after replace
%                                 testing, so code does not look for
%                                 '_TEST' file to replace)
%       'suffix',[char]         - manually add a suffix to the filename
%       'save_log'              - exports log of attempted saves with
%                                 success/failure status under
%                                 ~/CalcLogs/Saves_log_[starttime].txt
%       'version'               - set .mat file format to save data in. By
%                                 default, data is saved in v7.3 ([version]
%                                 = '-v7.3') which allows partial variable
%                                 loading but requires HDF5 loading
%                                 capability to import (which may cause
%                                 problems when attempting to open the file
%                                 in R or python). For more information,
%                                 see <a href =
%                                 "http://www.mathworks.com/help/matlab/ref/save.html?searchHighlight=v7.3&s_tid=doc_srchtitle#inputarg_fmt">here</a>
%                                 (official MathWorks summary) or <a href =
%                                 "http://undocumentedmatlab.com/blog/improving-save-performance">here</a>
%                                 (Undocumented Matlab blog post on the
%                                 pros and cons of using different file
%                                 formats (hint, v7 isn't usually useful,
%                                 but v6 may be)).
%
%       DEALING WITH DIFFERENT VARIABLE NAMING CONVENTIONS INSIDE NC FILES
%           'variable',[char]   - manually set the variable name to look
%                                 for in the source NetCDF file. By
%                                 default, the program looks for the
%                                 variable name set in Varnames.csv through
%                                 the input [VarIndices]. 
%           'lat_name',[char]   - manually set the name of the latitude
%                                 array. By default, the program looks for
%                                 'lat'. 
%           'lon_name',[char]   - manually set the name of the longitude
%                                 array. By default, the program looks for
%                                 'lon'. 
%           'lat_bnds_name',[char]  - manually set the name of the latitude
%                                     band array. By default, the program 
%                                     looks for 'lat_bnds'
%           'lon_bnds_name',[char]  - manually set the name of the
%                                     longitude band array. By default, the
%                                     program looks for 'lon_bnds'
%       NOTE: these above 5 flags are only supported with one variable
%       inputted (so length(VarIndices) == 1).
%
%   Saving convention:
%   [filevar]_[freq]_[model]_[exp]_[run]_[strtyr]_[endyr]_DATA.mat
%
%   NOTE: this function is part of the /project/moyer/ climate data file
%   ecosystem.
%
%   See also VARIOUS_DEFAULTS
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 07/25/2016

%% 1 Defaults
%Get defaults from various_defaults
various_defaults = matfile('various_defaults.mat');

%Saving conventions
years_to_save = 30; save_reg = true;
years_to_save_alt = 31; save_alt = true;
alt_folder = various_defaults.season_dir;
save_log = false;
skip_bounds = false;

replace_files = true;
filename_add = [];

%File variable names to search for 
var_search = [];
lat_search = []; %#ok<NASGU>
lon_search = []; %#ok<NASGU>
lat_bnds_search = []; %#ok<NASGU>
lon_bnds_search = []; %#ok<NASGU>

%Data conventions
version = '-v7.3';
expArray = {'rcp85','historical','piControl'};
dim4_val = [];
cust_file_strtyr = [];
cust_file_endyr = [];

%% 2 Misc setup
%Set regex delimiters to split up netcdf filenames
delim1 = '\_';
delim2 = '\.';
delim3= '\-';
delim = strcat(delim1,'|',delim2);

%% 3 Set behavior of optional function flags
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'replace'}
                replace_files = varargin{in_idx+1};
            case {'experiments'}
                expArray = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'dim4'}
                if ~isa(varargin{in_idx+1},'char')
                    dim4 = 'plev';
                    dim4_val = varargin{in_idx+1};
                else
                    dim4 = varargin{in_idx+1};
                    dim4_val = varargin{in_idx+2};
                end
            case {'testing'}
                filename_add = '_TEST';
            case {'num_years'}
                years_to_save = varargin{in_idx+1};
                if isa(varargin{in_idx+2},'numeric')
                    years_to_save_alt = varargin{in_idx+2};
                else
                    save_alt = false;
                end
            case {'save_case'}
                if strcmp(varargin{in_idx+1},'std')
                    save_alt = false;
                elseif strcmp(varargin{in_idx+1},'alt')
                    save_reg = false;
                elseif strcmp(varargin{in_idx+1},'all')
                else
                    warning('Saves:BadCaseInput',[varargin{in_idx+1},' is not a valid save case. Choose "std" or "alt"']);
                end
            case {'skip_bounds'}
                skip_bounds = varargin{in_idx+1};
            case {'save_log'}
                save_log = true;
            case {'year_range'}
                cust_file_strtyr = varargin{in_idx+1};
                cust_file_endyr = varargin{in_idx+2};
                save_alt = false;
                years_to_save = cust_file_endyr-cust_file_strtyr+1;
            case {'version'}
                version = varargin{in_idx+1};
            case {'suffix'}
                filename_add = varargin{in_idx+1};
            case {'variable'}
                var_search = varargin{in_idx+1};
            case {'lon_name'}
                lon_search = varargin{in_idx+1}; %#ok<NASGU>
            case {'lat_name'}
                lat_search = varargin{in_idx+1}; %#ok<NASGU>
            case {'lon_bnds_name'}
                lon_bnds_search = varargin{in_idx+1}; %#ok<NASGU> 
            case {'lat_bnds_name'}
                lat_bnds_search = varargin{in_idx+1}; %#ok<NASGU>
        end
    end
end

%Pre-allocate log array
if save_log
    complete_log = cell(length(VarIndices)*length(modelArray),length(expArray));
end
%% 4 Process
%Set variable
for i = 1:length(VarIndices)
    %Define file variable identifiers (filevar == netCDF raw variable name,
    %filevarFN == post-processing variable name. These are not necessarily
    %the same; notably when different elevation bands are included in the
    %netCDF file, but only one elevation band per file is processed into
    %_DATA.mat files)
    [filevar,filevarFN,freq,~,~,~,~,~] = var_chars(VarIndices(i));
    %Set model
    for j = 1:length(modelArray)
        %Set experiment
        for experiment = 1:length(expArray)
            %Save current variable state for clearing workspace purposes
            initial_vars = who;
            try
                %% 4.1 Identify netCDF files for defined variable
                %Get all netCDF files for defined variable
                FileNameSearchString=strcat(various_defaults.raw_data_dir,modelArray{j},'/',filevar,freq,modelArray{j},'_',expArray{experiment},'*.nc');
                FileNames = dir(FileNameSearchString);
                num_files = length(FileNames);
                if num_files==0
                    error('Saves:NoFiles',['No files matching ',FileNameSearchString,' found, no saves processed for ',filevarFN,freq,modelArray{j}])
                end
                
                %% 4.2 Get characteristics of variable and variable file structure
                %Get calendar type
                try
                    Calendar = ncreadatt([various_defaults.raw_data_dir,modelArray{j},'/',FileNames(1).name],'time','calendar');
                catch
                    prompt_str = ['The global attribute "calendar" was not found in ',...
                        [various_defaults.raw_data_dir,modelArray{j},'/',FileNames(1).name],...
                        '. Please manually select one of the following supported calendar formats (CMIP5 conventions) by inserting the corresponding integer:\n',...
                        '[1]: 365_day (choose this if the data is monthly or lower in frequency for the simplest function run) \n',...
                        '[2]: 360_day\n',...
                        '[3]: gregorian\n',...
                        '[4]: proleptic_gregorian\n',...
                        '[5]: no_leap\n',...
                        '[6]: julian\n'];
                    prompt_str = [prompt_str,'Insert a number to continue with the corresponding variable, or 0 to break.']; %#ok<AGROW>
                    calendars = {'365_day','360_day','gregorian','proleptic_gregorian','no_leap','julian'};
                    var_id = input(fprintf(prompt_str));
                    if var_id == 0
                        error('Save:varbreak',['Save run for ',filenameS,' interrupted.'])
                    else
                        Calendar = calendars{var_id};
                    end
                end
                %Pre-allocate Memory for Cell Arrays
                date_ranges=cell(1,num_files);
                run=cell(1,num_files);
                %Deconstruct netcdf filename into file characteristics
                for f=1:num_files
                    Varstr = FileNames(f).name;
                    splitStr = regexp(Varstr,delim,'split');
                    date_ranges{f} = splitStr{6};
                    run{f} = splitStr{5};
                end
                
                %% 4.3 Decide at which file to begin loading data
                %Get the last ending year of data for this var, model, exp
                %contained in filesystem (over arbitrary number of .nc
                %files), assuming date form of YYYY*, but also supporting
                %year-only dates of the form Y{1,4}.
                yrrSplit_end = regexp(date_ranges{num_files},delim3,'split');
                enddate = yrrSplit_end{2}; files_endyr = str2double(enddate(1:min(length(yrrSplit_end{2}),4)));
                
                %Get start years of each .nc file for this var, model, exp
                file_strtyrs = zeros(num_files,1);
                for f = 1:num_files
                    yrrSplit_strt = regexp(date_ranges{f},delim3,'split');
                    file_strtyrs(f) = str2double(yrrSplit_strt{1}(1:1:min(length(yrrSplit_strt{1}),4)));
                end
                %Get index of the file with the latest startdate that still
                %allows the full desired data length to be extracted
                [~,first_file] = max(file_strtyrs(files_endyr-file_strtyrs>max(years_to_save,years_to_save_alt)));
                if isempty(first_file)
                    first_file = 1;
                    warning('Saves:NotEnoughFiles',[num2str(max(years_to_save,years_to_save_alt)),' years of netcdf files not found. ',...
                        'The maximum found date range is ',num2str(file_strtyrs(1)),'-',num2str(files_endyr),...
                        '. This range will be used for all saves for ',filevarFN,freq,' ',modelArray{j},' ',expArray{experiment}])
                    years_to_save = files_endyr-file_strtyrs(1)+1;
                    years_to_save_alt = files_endyr-file_strtyrs(1)+1;
                end
                
                %Set files start date from index found above to minimize
                %imported files
                yrrSplit_strt = regexp(date_ranges{first_file},delim3,'split');
                strtdate = yrrSplit_strt{1}; files_strtyr = str2double(strtdate(1:min(length(yrrSplit_strt{1}),4)));
                
                %Truncate filenames to just those needed
                run = run(first_file:end);
                date_ranges = date_ranges(first_file:end);
                num_files = length(date_ranges); %updated number of (needed) files
                
                %Set last file startdate
                yrrSplit_strt = regexp(date_ranges{num_files},delim2,'split');
                strtdate = yrrSplit_strt{1}; last_strtyr = str2double(strtdate(1:min(length(yrrSplit_strt{1}),4)));
                
                %% 4.4 Get start / end indices for desired data length
                %Get number of data points in a year
                %if strfind(freq,'day')>0 
                if regexp(freq,'.[Dd]ay.')>0 %Daily data (julian/greg. are turned into 365_day below)
                    if strcmpi(Calendar,'noleap') || strcmpi(Calendar,'365_day') || strcmpi(Calendar,'gregorian') || strcmpi(Calendar,'proleptic_gregorian') || strcmpi(Calendar,'julian')
                        subyear_units = 365;
                    elseif strcmpi(Calendar,'360_day')
                        subyear_units = 360;
                    end
                elseif strfind(freq,'mon')>0 %Monthly data
                    subyear_units = 12;
                elseif contains(freq,'ann') || contains(freq,'yrr') %Yearly data
                    subyear_units = 1;
                else
                    error('Saves:NoFreqSupport',['Saves.m currently does not support files saved at frequency ',freq])
                end
                
                %Identify leap years (only for daily data)
                if strfind(freq,'day')>0
                    if strcmpi(Calendar,'noleap') || strcmpi(Calendar,'365_day') || strcmpi(Calendar,'360_day')
                        %(if no leap years, then the leap_year vector is just 0s)
                        leap_years = zeros(files_endyr-files_strtyr+1,1);
                    elseif strcmpi(Calendar,'julian')
                        years = files_strtyr:files_endyr;
                        %Identify leap years using standard julian leap year
                        %algorithm
                        leap_years = ~rem(years,4)~=0;
                    else
                        years = files_strtyr:files_endyr;
                        %Identify leap years using standard gregorian leap year
                        %algorithm
                        leap_years = ~(rem(years,4)~=0 | (rem(years,4)==0 & rem(years,100)==0 & rem(years,400)~=0));
                    end
                else
                    leap_years = zeros(files_endyr-files_strtyr+1,1);
                end
                
                %Process irregular start times
                if length(strtdate)<5
                    sub_startbias = 0;
                    if ~contains(freq,'yrr') && ~contains(freq,'ann') %Y{1-4} is fine for annual data, so no warning is necessary
                        warning('SAVES:ShortDate',['The start date of the time series for ',...
                            ' implied by the filename is less than 5 characters long, suggesting it is of the format Y{1-4}. ',...
                            'Accordingly, it is assumed the time series starts on January 1st / January / etc., and the "startbias" ',...
                            '(the number of timesteps skipped to correct for non-traditional start dates) is set to 0. Please verify if that feels fishy.'])
                    end
                elseif (contains(freq,'day') && strcmp(strtdate(5:8),'0101')) || (contains(freq,'mon') && strcmp(strtdate(5:6),'01'))
                    sub_startbias = 0;
                else
                    if strfind(freq,'day')>0 %Daily data
                        %If the file starts on a non-traditional (not
                        %01/01) day, start indexing on the next 01/01)
                        if strcmpi(Calendar,'360_day') %yearfrac is a financial toolbox function, might not always work
                            sub_startbias = yearfrac(strtdate(5:8),'1230',2)*subyear_units; %2 is for actual/360-day calendar
                        else
                            sub_startbias = yearfrac(strtdate(5:8),'1231',3)*subyear_units; %3 is for actual/365-day calendar
                            %If the first year is a leap year and the file
                            %starts at a non-traditional starting date
                            %before Feb 29th, add that day to the (to be
                            %skipped) start bias
                            if leap_years(1) && str2double(strtdate(5:8)) < 31 %representing 03/01, aka after any potential leap year)
                                sub_startbias = sub_startbias+1;
                            end
                        end
                    else %Monthly data
                        sub_startbias = subyear_units-num2str(strtdate(5:6))+1;
                    end
                    %Adjust leap_years array, since now effectively
                    %starting a year later (that leap day is dealt with
                    %below)
                    leap_years = leap_years(2:end);
                end
                
                %Process end year (year_endbias is the number of file years
                %to skip at the end of the last file)
                if ~isempty(cust_file_endyr)
                    year_endbias = cust_file_endyr-files_endyr;
                elseif files_endyr==2100 %By convention, save till 2099 (for consistency; some models only saved RCP runs to 2099)
                    year_endbias = -1;
                else
                    year_endbias = 0;
                end
                if year_endbias > 0 
                    error('Saves:FilesTooShort',['Though ',num2str(cust_file_endyr),' is desired as a final year for processed data, saved netCDF files only extend to ',num2str(files_endyr),'. Please pick an earlier end year or use default settings.'])
                end
                
                %Process start year (year_startbias is the number of file
                %years to skip before beginning import)
                if isempty(cust_file_strtyr)
                    year_startbias = ((files_endyr+year_endbias)-max(years_to_save,years_to_save_alt)+1)-files_strtyr;
                else %If custom start
                    year_startbias = cust_file_strtyr-files_strtyr-1;
                end
                if year_startbias < 0 
                    error('Saves:FilesTooShort',['Though ',num2str(cust_file_strtyr),' is desired as a first year for processed data, saved netCDF files only begin at ',num2str(files_strtyr),'. Please pick a later start year or use default settings.'])
                end
                
                %% 4.5 Define final (to save) filenames
                %Since year_startbias lists the number of complete years to
                %skip, the filename needs to take into account the first
                %full year that it contains (in other words,
                %year_startbias+1)
                if year_startbias~=0
                    std_filestrt = files_strtyr+year_startbias+1; 
                else 
                    std_filestrt = files_strtyr; 
                end
                
                filenameS = [various_defaults.raw_data_dir,modelArray{j},'/',filevarFN,freq,modelArray{j},'_',expArray{experiment},'_',run{1},'_',num2str(std_filestrt),'-',num2str(files_endyr+year_endbias),'_DATA.mat'];
                filename_alt = [various_defaults.proc_data_dir,modelArray{j},alt_folder,filevarFN,freq,modelArray{j},'_',expArray{experiment},'_',run{1},'_',num2str(std_filestrt+(years_to_save-years_to_save_alt)),'-',num2str(files_endyr+year_endbias),'_DATA.mat'];
                
                %% 4.6 Set which conventions to save
                %By default
                calc_reg = save_reg; calc_alt = save_alt;
                %If only calculate files that don't already exist,...
                if ~replace_files
                    if exist(filenameS,'file')
                        calc_reg = false;
                    end
                    if exist(filename_alt,'file')
                        calc_alt = false;
                    end
                end
                %Only run process if at least one of the two files defined
                %by the above filenames is needed
                if ~calc_reg && ~calc_alt
                    calc = false;
                else
                    calc = true;
                end
                
                %% 4.7 Calculate, if needed
                if calc
                    %% 4.7.1 Define dimension counts and indexes
                    file_info = ncinfo([various_defaults.raw_data_dir,modelArray{j},'/',filevar,freq,modelArray{j},'_',expArray{experiment},'_',run{1},'_',date_ranges{1},'.nc']);
                    %Make sure the filevar exists in the nc files
                    %(interior_var is the variable name on the "interior"
                    %of the file, as opposed to filevar, which is the
                    %variable name in the filename. These shouldn't be
                    %different, but occasionally are). 
                    if isempty(var_search)
                        if ~isempty(structfind(file_info.Variables,'Name',filevar))
                            interior_var = filevar;
                        else
                            prompt_str = ['The variable "',filevar,'" was not found in ',...
                                various_defaults.raw_data_dir,modelArray{j},'/',FileNames(1).name,...
                                '. Please manually select one of the following variables to load by inserting the corresponding integer:\n'];
                            for nc_var_idxs = 1:length(file_info.Variables)
                                prompt_str = [prompt_str,'[',num2str(nc_var_idxs),']: ',file_info.Variables(nc_var_idxs).Name,'\n']; %#ok<AGROW>
                            end
                            prompt_str = [prompt_str,'Insert a number to continue with the corresponding variable, or 0 to break. ',...
                                'If this program is being run with many variables / models, I suggest breaking the code run and ',...
                                're-running with the flag "variable" set to the variable name actually present in the file to avoid ',...
                                'having to sit by this program while it is interrupted every time it runs into this model ',...
                                '(note: this flag only works with one variable and set of models with the same unconventional variable naming scheme at a time).']; %#ok<AGROW>
                            var_id = input(fprintf(prompt_str));
                            if var_id == 0
                                error('Save:varbreak',['Save run for ',filenameS,' interrupted.'])
                            else
                                interior_var = file_info.Variables(var_id).Name;
                            end
                        end
                    else
                        %If specifying a file-interior variable name, set
                        %program to search for it. 
                        interior_var = var_search;
                    end
                    
                    %Get number of dimensions (lat, lon, time, etc.) for var
                    var_size = file_info.Variables(structfind(file_info.Variables,'Name',interior_var)).Size;
                    
                    %Get which dimension of the variable in the netcdf is time
                    time_dim = structfind(file_info.Variables(structfind(file_info.Variables,'Name',interior_var)).Dimensions,'Name','time');
                    
                    %Get base count starts and ends (ones = from beginning of
                    %that dimension in file, inf = till end of that dim)
                    %(are row vectors to comply with ncread format)
                    idx_strt1 = ones(1,length(var_size));
                    idx_strtf = ones(1,length(var_size));
                    idx_count1 = inf(1,length(var_size));
                    idx_countf = inf(1,length(var_size));
                    idx_interm_i = ones(1,length(var_size));
                    idx_interm_count = inf(1,length(var_size));
                    
                    %Get which dimension of the variable in the netcdf is
                    %elevation/other desired fourth dimension variable and
                    %insert it into count var
                    if ~isempty(dim4_val)
                        dim4_idx = structfind(file_info.Variables(structfind(file_info.Variables,'Name',interior_var)).Dimensions,'Name',dim4);
                        idx_strt1(dim4_idx) = dim4_val;
                        idx_strtf(dim4_idx) = dim4_val;
                        idx_count1(dim4_idx) = 1;
                        idx_countf(dim4_idx) = 1;
                        idx_interm_i(dim4_idx) = dim4_val;
                        idx_interm_count(dim4_idx) = 1;
                        var_size(dim4_idx) = 1; %Make sure indices treat the fourth dimension as singleton
                    end
                    
                    %Get index to start on on first file
                    idx_strt1(time_dim) = sub_startbias+year_startbias*subyear_units+sum(leap_years(1:year_startbias))+1;
                    %Get index to end on on last file
                    idx_countf(time_dim) = (files_endyr-last_strtyr+1+year_endbias)*subyear_units+sum(leap_years(length(leap_years)-(files_endyr-last_strtyr+1+year_endbias)+1:(length(leap_years)+year_endbias)));
                    %If only one netcdf file, subtract idx_strt1 since
                    %ncread works with counts, not indices
                    if num_files == 1
                        idx_countf(time_dim) = idx_countf(time_dim) - idx_strt1(time_dim)+1;
                        idx_count1 = idx_countf;
                    end
                    
                    %% 4.7.2 Load variable data
                    %Pre-allocate cells
                    filename = cell(num_files,1);
                    Raw = cell(num_files,1);
                    
                    %Load first file, starting at startbias
                    filename{1} = [various_defaults.raw_data_dir,modelArray{j},'/',filevar,freq,modelArray{j},'_',expArray{experiment},'_',run{1},'_',date_ranges{1},'.nc'];
                    Raw{1} = ncread(filename{1},filevar,idx_strt1,idx_count1);
                    %Fully load intermediate files
                    for f = 2:num_files-1
                        filename{f} = [various_defaults.raw_data_dir,modelArray{j},'/',filevar,freq,modelArray{j},'_',expArray{experiment},'_',run{f},'_',date_ranges{f},'.nc'];
                        Raw{f} = ncread(filename{f},filevar,idx_interm_i,idx_interm_count);
                    end
                    %Load last file, ending at endbias
                    if num_files ~= 1
                        filename{num_files} = [various_defaults.raw_data_dir,modelArray{j},'/',filevar,freq,modelArray{j},'_',expArray{experiment},'_',run{num_files},'_',date_ranges{num_files},'.nc'];
                        Raw{num_files} = ncread(filename{num_files},filevar,idx_strtf,idx_countf);
                    end
                    
                    %Convert into matrix
                    Raw = cat(time_dim,Raw{:});
                    %Add Raw to output struct
                  
                    
                    %% 4.7.3 Process Leap Years, if necessary (by removing leap-added day from end of that year)
                    idxs_tmp = cell((files_endyr-files_strtyr-year_startbias),1);
                    %Seed with first year (minus leap-added day, if extant)
                    idxs_tmp{1} = 1:subyear_units;
                    %For subsequent years, start index at next value from last
                    %of previous in common years, last value + 1 for the
                    %preceding being a leap year (thereby skipping it)
                    for y = 2:(files_endyr-files_strtyr+1-year_startbias+year_endbias)
                        idxs_tmp{y} = (max(idxs_tmp{y-1})+leap_years(year_startbias+y-1)+1):(max(idxs_tmp{y-1})+leap_years(year_startbias+y-1))+subyear_units;
                    end
                    %Turn indices into single vector
                    idxs_tmp = [idxs_tmp{:}];
                    %Get indices of all other dimensions (so say 1:nlon,
                    %1:nlat, etc. to keep those counts constant while taking
                    %out some time dimension points)
                    idxs = cell(length(var_size),1);
                    for dim = 1:length(var_size)
                        idxs{dim} = 1:var_size(dim);
                    end
                    idxs{time_dim} = idxs_tmp;
                    
                    %Clip out leap-added days
                    Raw = Raw(idxs{:});
                    
                    %Remove singleton dimensions (i.e. set 4th dim
                    %elevation)
                    Raw = squeeze(Raw);
                    
                    %% 4.7.4 Load Geographic Variables 
                    %Load latitude, longitude, latitude bounds, and
                    %longitude bounds, with support for files in which
                    %these are not called lat, lon, lat_bnds, and lon_bnds,
                    %respectively
                    if skip_bounds
                        aux_vars = {'lat','lon'};
                    else
                        aux_vars = {'lat','lon','lat_bnds','lon_bnds'};
                    end
                    for aux_var_idx = 1:length(aux_vars)
                        if isempty(var_search)
                            if ~isempty(structfind(file_info.Variables,'Name',aux_vars{aux_var_idx}))
                                outputs.(aux_vars{aux_var_idx}) = ncread(filename{1},aux_vars{aux_var_idx});
                            else
                                prompt_str = ['The variable "',aux_vars{aux_var_idx},'" was not found in ',...
                                    various_defaults.raw_data_dir,modelArray{j},'/',FileNames(1).name,...
                                    '. Please manually select one of the following variables to load as "',aux_vars{aux_var_idx},'" by inserting the corresponding integer:\n'];
                                for nc_var_idxs = 1:length(file_info.Variables)
                                    prompt_str = [prompt_str,'[',num2str(nc_var_idxs),']: ',file_info.Variables(nc_var_idxs).Name,'\n']; %#ok<AGROW>
                                end
                                prompt_str = [prompt_str,'Insert a number to continue with the corresponding variable, or 0 to break. ',...
                                    'If this program is being run with many variables / models, I suggest breaking the code run and ',...
                                    're-running with the flag "',aux_vars{aux_var_idx},'_name" set to the variable name actually present in the file to avoid ',...
                                    'having to sit by this program while it is interrupted every time it runs into this model ',...
                                    '(note: this flag only works with one variable and set of models with the same unconventional variable naming scheme at a time).']; %#ok<AGROW>
                                var_id = input(fprintf(prompt_str));
                                if var_id == 0
                                    error('Save:varbreak',['Save run for ',filenameS,' interrupted.'])
                                else
                                    outputs.(aux_vars{aux_var_idx}) = ncread(filename{1},file_info.Variables(var_id).Name);
                                end
                            end
                        else
                            %If specifying a file-interior variable name, set
                            %program to search for it.
                            outputs.(aux_vars{aux_var_idx}) = ncread(filename{1},eval([aux_vars{aux_var_idx},'_search']));
                        end
                    end
                    
                    %% 4.7.5 Output struct?
                    outputs.Raw = Raw; clear Raw
                    
                    %% 4.7.6 Save files
                    %Save longer time series
                    if years_to_save_alt >= years_to_save && calc_alt
                        if ~exist([various_defaults.proc_data_dir,modelArray{j},alt_folder],'file')
                            mkdir([various_defaults.proc_data_dir,modelArray{j},alt_folder])
                            disp([various_defaults.proc_data_dir,modelArray{j},alt_folder,' did not yet exist, and was created.'])
                        end
                        save([filename_alt,filename_add],'-struct','outputs',version)
                        disp([filename_alt,filename_add,' ',version,' saved'])
                        calc_alt = false; %make sure not saved twice
                    elseif calc_reg
                        save([filenameS,filename_add],'-struct','outputs',version)
                        disp([filenameS,filename_add,' ',version,' saved'])
                        calc_reg = false; %make sure not saved twice
                    end
                    
                    %Prune idxs down to smaller of two save conditions
                    % THIS SHOULD BE CHANGED, CURRENTLY ASSUMES RAW HAS 3
                    % DIMENSIONS< WHICH IS NOT ALWAYS THE CASE
                    idxs_tmp = abs(years_to_save_alt-years_to_save)*subyear_units+1:size(outputs.Raw,3);
                    idxs = cell(ndims(outputs.Raw),1);
                    for dim = 1:length(outputs.Raw)
                        idxs{dim} = 1:size(outputs.Raw,dim);
                    end
                    idxs{3} = idxs_tmp;
                    outputs.Raw = outputs.Raw(idxs{:});

                    
                    %Save shorter time series
                    if years_to_save >= years_to_save_alt && calc_alt
                        if ~exist([various_defaults.proc_data_dir,modelArray{j},alt_folder],'file')
                            mkdir([various_defaults.proc_data_dir,modelArray{j},alt_folder])
                            disp([various_defaults.proc_data_dir,modelArray{j},alt_folder,' did not yet exist, and was created.'])
                        end
                        save([filename_alt,filename_add],'-struct','outputs',version)
                        disp([filename_alt,filename_add,' ',version,' saved'])
                    elseif calc_reg
                        save([filenameS,filename_add],'-struct','outputs',version)
                        disp([filenameS,filename_add,' ',version,' saved'])
                    end
                    
                    %Store success message in log
                    if save_log
                        save_log_msg = [filevarFN,freq,modelArray{j},expArray{experiment},' Complete'];
                        complete_log{(i-1)*length(modelArray)+j,experiment} = save_log_msg;
                    end
                    
                else
                    disp(['File for ',modelArray{j},' ',filevarFN,freq,expArray{experiment},' already exists, was not replaced']);
                    %Store message that no figure needed in log
                    if save_log
                        save_log_msg = [filevarFN,freq,modelArray{j},expArray{experiment},' already exists'];
                        complete_log{(i-1)*length(modelArray)+j,experiment} = save_log_msg;
                    end
                end
                clearvars('-except',initial_vars{:});
                
            catch ME
                disp(ME.message);
                disp(ME.stack);
                %Store warning in log
                if save_log
                    save_log_wrn = ['****',filevarFN,freq,modelArray{j},expArray{experiment},' Not Complete****'];
                    complete_log{(i-1)*length(modelArray)+j,experiment} = save_log_wrn;
                end
                warning('Saves:Fail',['Filesaves for ',modelArray{j},' ',expArray{experiment},' ',filevarFN,freq,' Failed'])
                clearvars('-except',initial_vars{:});
            end
        end
        disp(['All files for ',modelArray{j},' ',filevarFN,' processed.']);
    end
end

%% 5 Export log as table, if desired
if save_log
    complete_log = complete_log';
    complete_log = complete_log(:);
    export_log = cell2table(complete_log','VariableNames',{'Saves_attempted_executions'});
    writetable(export_log,['/home/kschwarzwald/CalcLogs/Saves_log_',num2str(startTimestamp(1)),'-',num2str(startTimestamp(2)),'-',num2str(startTimestamp(3)),'-',num2str(startTimestamp(4)),'-',num2str(startTimestamp(5)),'-',num2str(startTimestamp(6)),'.txt']);
end
