function varargout = name_chars(model,expArray,filevar,freq,varargin)
% NAME_CHARS    get filename characteristics
%   [yrr,run,strtyr,endyr] = NAME_CHARS(model,expArray,filevar,freq)
%   extracts the year range, experiment run, and file start and end dates
%   from the filename of the first (alphabetically listed) _DATA.mat file
%   in the directory defined by [model], the experiment(s) in [expArray],
%   the variable identified [filevar], and the frequency [freq]. All
%   outputs are length(expArray) x 1 cell arrays, in the order given by
%   expArray. [yrr] gives the full [dates] segment of the filename in a
%   string, [run] a string of the file run, and [strtyr] and [endyr]
%   provide (numeric) starting and ending years of the file. If two
%   experiments are listed in [expArray], NAME_CHARS will attempt to match
%   the length of the two experiment runs by cycling through permutations
%   of saved run lengths.
%
%   If more than one experiment desired, [expArray] is a num_exps x 1 cell
%   array, else it can be inputted as a string. All other inputs are
%   strings.
%
%   [yrr,run,strtyr,endyr,num_conventions] =
%   NAME_CHARS(model,expArray,filevar,freq) also outputs the number of
%   files that fit the search given below, in a numerical length(expArray)
%   x 1 array.
%
%   By default, NAME_CHARS searches for a _DATA.mat file in the directory
%   [raw_data_dir]/[model]/ of the form
%   "[var]_[freq]_[model]_[experiment]_[run]_[dates]_DATA.mat".
%   or in [proc_data_dir]/[model]/[season_dir] for files of the
%   form
%   "[var]_[freq]_[model]_[experiment]_[run]_[dates]_sqrtPower_[seas].mat"
%   for seasonal data.
%
%   NAME_CHARS(...,'[flag]',[params],...) modify program run as below:
%       'directory',[dir]    - searches the directory given by the string
%                              [dir] rather than the default location.
%       'use_convention',[#] - returns characteristics of the #th
%                              (alphabetically) file. Default is 1. If
%                              use_convention > num_conventions, outputs
%                              are empty and a warning is thrown.
%       'season',[season],[season_folder]
%                            - adds '*[season]' to the file search string.
%                              Optionally change the season folder (in
%                             [proc_data_dir][model]) from which to
%                              get and save files by also setting
%                              [season_folder] = '/[folder name]/' (by
%                              default [season_dir]). Setting season to
%                              'all' does not change the running of the
%                              program.
%       'start_year',[num/cell] - searches for a specific start year
%                                 (specifically from the first portion the
%                                 next part of the CMIP5 standard filename
%                                 after the [run] slot). By default is set
%                                 to [] (not searched for). If more than
%                                 one set of file characteristics is
%                                 outputted, this can be a cell array (in
%                                 which case each element corresponds to
%                                 searching for each set of file
%                                 characteristics); otherwise it is assumed
%                                 the desired start year is the same for
%                                 both.
%       'end_year',[num/cell]   - searches for a specific end year
%                                 (specifically from the last portion the
%                                 next part of the CMIP5 standard filename
%                                 after the [run] slot). By default is set
%                                 to [] (not searched for). If more than
%                                 one set of file characteristics is
%                                 outputted, this can be a cell array (in
%                                 which case each element corresponds to
%                                 searching for each set of file
%                                 characteristics); otherwise it is assumed
%                                 the desired end year is the same for
%                                 both.
%       'run',[str/cell]        - also searches for a specific run
%                                 (specifically the next part of the CMIP5
%                                 standard filename after the [experiment]
%                                 slot. By default is set to [], an empty
%                                 string, and is therefore not explicitly
%                                 searched for). If more than
%                                 one set of file characteristics is
%                                 outputted, this can be a cell array (in
%                                 which case each element corresponds to
%                                 searching for each set of file
%                                 characteristics); otherwise it is assumed
%                                 the desired run is the same for both.
%       'file_end',[str]     - searches for files that end in '[str].mat'
%                              instead of default '*DATA.mat'.
%       'match_length',[log] - set whether to attempt to match length of
%                              the two experiment runs (def: true if >= 2
%                              experiment runs are inputted. If set to true
%                              and more than two experiment runs are
%                              inputted, code will attempt to jointly match
%                              all three, and output highest convention (as
%                              returned in alphanumeric order by system)
%                              indexed combination of all experiments if no
%                              match is found).
%
%   Example:
%   [___] = name_chars('NorESM1-M',{'rcp85','piControl'},'ta850','day')
%   will return:
%       yrr = {'2071-2099','830-858'}
%       run = {'r1i1p1_','r1i1p1_'}
%       strtyr = {2071,830};
%       endyr = {2099,858};
%       num_conventions = {2,1};
%   (A warning will be thrown saying multiple filenames were found for
%   ta850 rcp85, as can be seen in num_conventions{1}>1. Note also that in
%   deciding between the two conventions, characteristics for the one
%   matching total file length to the piControl run are outputted)
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also VAR_CHARS
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 07/21/2017

various_defaults = matfile('various_defaults.mat');
diff_directory = false;
use_convention = ones(1,length(expArray));
file_end = '*DATA.mat'; cust_file_end = false;
season = 'all'; season_folder = various_defaults.season_dir;
run_set_tmp = [];
strtyr_set_tmp = [];
endyr_set_tmp = [];


%If two experiments, by default attempt to match lengths
if length(expArray)>=2
    match_length = true;
else
    match_length = false;
end

%If expArray is inputted as a single string, convert to cell
if isa(expArray,'char')
    expArray = {expArray};
end

%If freq does not have underscores, add them
if isempty(strfind(freq,'_'))
    freq = ['_',freq,'_'];
end

%Set behavior of optional function flags
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'directory'}
                diff_directory = true;
                alt_directory = varargin{in_idx+1};
            case {'use_convention'}
                for experiment = 1:length(expArray)
                    use_convention(experiment) = varargin{in_idx+experiment};
                end
            case {'season'}
                season = varargin{in_idx+1};
                if ~strcmp(season,'all')
                    file_end = ['*',season,'.mat'];
                    if length(varargin)>in_idx+1
                        if strcmp(varargin{in_idx+2}(1),'/')
                            season_folder = varargin{in_idx+2};
                        end
                    end
                end
            case {'file_end'}
                file_end = varargin{in_idx+1}; cust_file_end = true;
            case {'match_length'}
                match_length = varargin{in_idx+1};
            case {'start_year'} 
                strtyr_set_tmp = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'end_year'} 
                endyr_set_tmp = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'run'}
                run_set_tmp = varargin{in_idx+1}; varargin{in_idx+1} = 0;
                
        end
    end
end
%Set regex delimiters to split up netcdf filenames
delim1 = '\_';
delim2 = '\.';
delim3= '\-';
delim = strcat(delim1,'|',delim2);

%Set directory of searching by season
if strcmp(season,'all')
    seas_dir = '/';
    main_directory = various_defaults.raw_data_dir;
else
    seas_dir = season_folder;
    main_directory = various_defaults.proc_data_dir;
    if ~cust_file_end
        file_end = ['*sqrtPower_',season,'.mat'];
    end
end

%Set final directory to look in
if ~diff_directory
    directory = strcat(main_directory,model,seas_dir);
else
    directory = alt_directory;
end

%Preallocate arrays
num_exps = length(expArray);
yrr = cell(1,num_exps);
run = cell(1,num_exps);
strtyr = cell(1,num_exps);
endyr = cell(1,num_exps);
num_conventions = ones(1,num_exps);
for experiment = 1:num_exps
    %Deal with both cell and individual inputs
    if isa(strtyr_set_tmp,'cell')
        strtyr_set = strtyr_set_tmp{experiment};
    else
        strtyr_set = strtyr_set_tmp;
    end
    if isa(endyr_set_tmp,'cell')
        endyr_set = endyr_set_tmp{experiment};
    else
        endyr_set = endyr_set_tmp;
    end
    if isa(run_set_tmp,'cell')
        %Add underscore to the run name to keep filename stable
        if ~isempty(run_set_tmp{experiment})
            run_set = [run_set_tmp{experiment},'_'];
        else
            run_set = [];
        end
    else
        if ~isempty(run_set_tmp)
            run_set = [run_set_tmp,'_']; 
        else
            run_set = [];
        end
    end
    
    %Set year string if desired
    if isempty(strtyr_set) && isempty(endyr_set)
        year_set = [];
    else
        if isempty(strtyr_set) && ~isempty(endyr_set)
            year_set = ['*-',num2str(endyr_set)];
        elseif ~isempty(strtyr_set) && isempty(endyr_set)
            year_set = ['*',num2str(strtyr_set),'-*'];
        else
            year_set = ['*',num2str(strtyr_set),'-*',endyr_set];
        end
    end
    
    FileNameSearchString = [directory,filevar,freq,model,'_',expArray{experiment},'_',run_set,year_set,file_end];
    FileNames = dir(FileNameSearchString);
    num_conventions(experiment) = length(FileNames);
    if (use_convention(experiment) > num_conventions(experiment)) || use_convention(experiment)==0
        NoExistMsg = ['use_convention is greater than number of conventions, ',...
            FileNameSearchString,' does not exist'];
        warning('name_chars:FileNotFound',NoExistMsg)
        yrr{experiment} = '';
        run{experiment} = '';
        strtyr{experiment} = '';
        endyr{experiment} = '';
    else
        
        %Get run and timeframe characteristics for file
        Varstr = FileNames(use_convention(experiment)).name;
        splitStr = regexp(Varstr,delim,'split');
        yrr{experiment} = splitStr{6};
        run{experiment} = strcat(splitStr{5},'_');
        
        %Define start and end years for files
        yrrSplit = regexp(yrr,delim3,'split');
        strtdate = yrrSplit{experiment}(1);
        strtyr{experiment} = str2double(strtdate);
        enddate = yrrSplit{experiment}(2);
        endyr{experiment} = str2double(enddate);
        
        if num_conventions(experiment) > 1 && ~match_length
            warn_Msg = ['More than one ',file_end,' file for ',model,' ',expArray{experiment},...
                ' ',freq,' ',filevar,' ',run_set,' found, ',num2str(strtyr{experiment}),'-',num2str(endyr{experiment}),' name characteristics outputted for ',expArray{experiment}];
            warning('name_chars:MultConvs',warn_Msg)
        end
    end
end
if match_length
    if length(unique([endyr{:}]-[strtyr{:}]))>1
        if any(num_conventions>1); convs = cell(length(expArray),1);
            for experiment = 1:length(expArray)
                convs{experiment} = 1:num_conventions(experiment);
            end
            possible_convs = combvec(convs{:})';
            n = 1;
            while length(unique([endyr{:}]-[strtyr{:}]))>1
                %Break if reached past number of conventions
                if n > size(possible_convs,1)
                    warning('name_chars:LengthMismatch',...
                        [model,' ',strjoin(expArray,', '),' runs do not jointly have the same length. Highest convention index combination (',...
                        strjoin(cellfun(@(x) num2str(x),strtyr,'UniformOutput',0),', '),' to ',strjoin(cellfun(@(x) num2str(x),endyr,'UniformOutput',0),', '),...
                        ' for ',strjoin(expArray,', '),', respectively) outputted, further verification suggested.'])
                    break
                end
                conv_try = possible_convs(n,:); %Each row of permutations
                for experiment = 1:length(expArray) %Redo processing above
                    FileNameSearchString = [directory,filevar,freq,model,'_',expArray{experiment},'_',run_set,year_set,file_end];
                    FileNames = dir(FileNameSearchString);
                    %Get run and timeframe characteristics for file
                    Varstr = FileNames(conv_try(experiment)).name;
                    splitStr = regexp(Varstr,delim,'split');
                    yrr{experiment} = splitStr{6};
                    run{experiment} = strcat(splitStr{5},'_');
                    
                    %Define start and end years for files
                    yrrSplit = regexp(yrr,delim3,'split');
                    strtdate = yrrSplit{experiment}(1);
                    strtyr{experiment} = str2double(strtdate);
                    enddate = yrrSplit{experiment}(2);
                    endyr{experiment} = str2double(enddate);
                end
                n = n+1; 
            end
        else
            warning('name_chars:LengthMismatch',...
                        [model,' ',strjoin(expArray,', '),' runs do not jointly have the same length. Highest convention index combination (',...
                        strjoin(cellfun(@(x) num2str(x),strtyr,'UniformOutput',0),', '),' to ',strjoin(cellfun(@(x) num2str(x),endyr,'UniformOutput',0),', '),...
                        ' for ',strjoin(expArray,', '),', respectively) outputted, further verification suggested.'])
        end
    end
end

varargout{1} = yrr;
varargout{2} = run;
varargout{3} = strtyr;
varargout{4} = endyr;
if nargout == 5
    varargout{5} = num_conventions;
end


end