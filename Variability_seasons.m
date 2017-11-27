% VARIABILITY_SEASONS   calculate variability within and between seasons.
%
%   VARIABILITY_SEASONS(VarIndices,modelArray) calculates the contribution
%   of signals in various frequency bands to the standard deviation of
%   variables from models given by [VarIndices] and [modelArray], separated
%   by season. All data is first linearly detrended and deseasonalized (by
%   a harmonic filter for daily data and by simply removing the average of
%   all of each month from each of that month for monthly data). Then a
%   power spectrum is calculated for each individual season using a
%   standard FFT procedure, which is then averaged to produce a seasonal
%   power spectrum over the full timeframe of the data file.
%   Intraseasonal/intraannual standard deviations are calculated by taking
%   the square root of sums of the spectral density up to set frequencies,
%   scaled by the length of the season. Interseasonal/interannual
%   variability is calculated by taking the standard deviation of the
%   average of each detrended season. The local mean of each variable (for
%   the given season) is also stored.
%
%   The program outputs two files for each independent _DATA file and
%   season processed - one file _sqrtPower_[season] with a structure array
%   containing all variability data separated by frequency band, a word
%   description of each frequency band, the frequency band endpoints (in
%   base frequency units of the variable - days for daily or months for
%   monthly data), and an identifier (i.e. "XF" for interannual data), in
%   addition to (separately) the lat and lon variables from the original
%   file - and one file _LocalMeans_[season] with the local mean of the
%   variable (for the given season).
%
%   Default frequency bands (by default, 30 seasons or largest possible
%   subset taken from data):
%       ID      Daily       Monthly     Notes
%       HF      1-5         -           saved as [0,5] by convention
%       MF      3-15        -
%       LF      15-90       1-3
%       XF      > 1 yr      > 1 yr      saved as [12,360] or [365,36500] by
%                                       convention
%
%   Sample Output for 2-band monthly [_sqrtPower] file:
%   StdDevs.Name = {'< 3 months';'> 1 year'};
%   StdDevs.ID = {'LF';'XF'};
%   StdDevs.Size = {0,3;12,360};
%   StdDevs.Data = {[nlon x nlat array];[nlon x nlat array]};
%   lat = [nlat x 1 array]
%   lon = [nlon x 1 array]
%
%   Data requirements: [_DATA.mat] files for each variable, experiment,
%                      model in the folder given by [folder_seasons] (by
%                      default [raw_data_dir]/[model]/Seasons_DS6/)
%   Function requirements (on path): name_chars, var_chars
%
%   VARIABILITY_SEASONS(...,'[flag]',[params],...) modify program run as
%   below:
%
%       'season',[name](,...) - manually set season(s) to process (default
%                             {'JJA','DJF'}). If an array other than
%                             {'DJF'} or {'JJA'} is inputted (or if {'DJF'}
%                             or {'JJA'} are desired, but with non-default
%                             other season parameters), the season [name]
%                             should be followed by the following cell
%                             arrays for *BOTH* daily and monthly data,
%                             even if variables with only one frequency are
%                             chosen for processing (the following example
%                             is for processing 30 90-day/3-month seasons
%                             starting on February 1st. A possible [name]
%                             from above would be {'FMA'} for
%                             Feb-March-April):
%                   [length],        length of season (i.e. {90;3})
%                   [start],         season start date index (i.e. {32;2})
%                   [timeframe],     length in years of _DATA file wanted /
%                                    number of *calendar* years desired
%                                    (i.e. {30}) - if season does not cross
%                                    calendar years, this is equal to number
%                                    of seasons wanted. If it does, this is
%                                    equal to 1+number of seasons wanted.
%                                    If _DATA file has fewer than
%                                    [timeframe] years of data, the maximum
%                                    possible number of seasons is taken.
%                             Multiple seasons can be added simultaneously,
%                             with each season having its own column in the
%                             inputs above.
%
%       'experiments',[cell]- manually set experiments to process (default
%                             {'rcp85','historical','piControl'})
%       'num_harmonics',[#] - change the number of harmonics removed in
%                             deseasonalization process for daily data
%       'monthly_bands',[names],[lims],[ID]  - set new endpoints for
%                                         monthly data, default is [names]
%                                         = {'< 3 months'}, [lims] = [0,3],
%                                         and [ID] = {'LF','XF'}
%       'daily_bands',[names],[lims],[ID]    - set new endpoints for daily
%                                         data, default is [names] = {'< 3
%                                         days','3 - 15 days','15 - 90
%                                         days'}, [lims] =
%                                         [0,3;3,15;15,90], and [ID] =
%                                         {'HF','MF','LF'};
%       'multi_processor',[log],[#]     - set whether to use multiple
%                                         processors and set number of [#]
%                                         processors to be used (def: true,
%                                         12)
%       'batch_processing',[log]        - sets job storage location
%                                         explicitly at
%                                         /tmp/kschwarzwald/[SLURMJOBID]
%                                         (def: false; only useful in batch
%                                         processing)
%       'replace',[log]     - set whether to replace existing files (def:
%                             true)
%       'save_log'          - exports log of attempted saves with
%                             success/failure status under
%                        ~/CalcLogs/Variability_seasons_log_[starttime].txt
%
%
%   Saving convention:
%   [filevar]_[freq]_[model]_[exp]_[run]_[strtyr]_[endyr]_sqrtPower_[season].mat
%   [filevar]_[freq]_[model]_[exp]_[run]_[strtyr]_[endyr]_LocalMeans_[season].mat
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also VARIABILITY, VARIABILITY_CI
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 11/26/2017

function Variability_seasons(VarIndices,modelArray,varargin)

%Set clock (for logging purposes)
format shortG
startTimestamp = clock;
tic

%% Set default values and setup program
various_defaults = matfile('various_defaults.mat');
replace_files = true;
num_models = length(modelArray);
num_harmonics=12;
expArray = {'rcp85','historical','piControl'};
num_exps = length(expArray);
multi_proc = true; num_proc = 12; batch_processing = false;
%pc_jobstorloc = strcat('/tmp/kschwarzwald/', getenv('SLURM_JOB_ID'));
pc_jobstorloc = strcat(various_defaults.tmp_dir,randi(10000));
filename_add = [];
save_log = false;
folder = various_defaults.season_dir;

%Seasonal Setup
seasonArray = {'JJA','DJF'};
seasonLength = {90,90;3,3}; %Length (in days) of season
seasonStart = {153,336;6,12}; %Start day (of year) of season
prefTimeframe = {30,31}; %How many years of data preferred (DJF intended to be calculated using 31 years of data, since starting in december therefore gives 30 'winters' of data)

%Daily Freq Bands
daily_freq_bands = [0,3;3,15;15,90]; daily_num_bands = size(daily_freq_bands,1);
daily_freq_band_names = {'< 3 days','3 - 15 days','15 - 90 days'};
daily_freq_band_ID = {'HF','MF','LF'}; freq_idx_day = 1:daily_num_bands;

%Monthly Freq Bands
monthly_freq_bands = [0,3]; monthly_num_bands = size(monthly_freq_bands,1);
monthly_freq_band_names = {'< 3 months'};
monthly_freq_band_ID = {'LF'}; freq_idx_month = 1:monthly_num_bands;

%% Set behavior of optional function flags
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'season'}
                seasonArray_tmp = varargin{in_idx+1}; varargin{in_idx+1} = 0;
                seas_test_idx = find(strcmp(seasonArray,seasonArray_tmp{1}));
                if length(seasonArray_tmp) == 1 && ~isempty(seas_test_idx) && (length(varargin) == in_idx+1 || ~isa(varargin{in_idx+2},'cell')) %If desired season is DJF or JJA with default settings
                    seasonLength = seasonLength(:,seas_test_idx);
                    seasonStart = seasonStart(:,seas_test_idx);
                    prefTimeframe = prefTimeframe(seas_test_idx);
                else %If desired season not DJF or JJA, or is, but with different number of seasons, etc.
                    seasonArray = seasonArray_tmp; clear seasonArray_tmp
                    seasonLength = varargin{in_idx+2}; varargin{in_idx+2} = 0;
                    seasonStart = varargin{in_idx+3}; varargin{in_idx+3} = 0;
                    prefTimeframe = varargin{in_idx+4}; varargin{in_idx+4} = 0;
                end
            case {'replace'}
                replace_files = varargin{in_idx+1};
            case {'experiments'}
                expArray = varargin{in_idx+1}; varargin{in_idx+1} = 0;
                num_exps = length(expArray);
            case {'num_harmonics'}
                num_harmonics = varargin{in_idx+1};
                if isa(varargin{in_idx}+1,'integer') == 0 || varargin{in_idx+1}<0
                    error('num_harmonics must be a positive integer')
                end
            case {'monthly_bands'}
                monthly_freq_band_names = varargin{in_idx+1}; varargin{in_idx+1} = 0;
                monthly_freq_bands = varargin{in_idx+2}; varargin{in_idx+2} = 0;
                monthly_freq_band_ID = varargin{in_idx+3}; varargin{in_idx+3} = 0;
                monthly_num_bands = size(monthly_freq_bands,1); freq_idx_month = 1:monthly_num_bands;
            case {'daily_bands'}
                daily_freq_band_names = varargin{in_idx+1}; varargin{in_idx+1} = 0;
                daily_freq_bands = varargin{in_idx+2}; varargin{in_idx+2} = 0;
                daily_freq_band_ID = varargin{in_idx+3}; varargin{in_idx+3} = 0;
                daily_num_bands = size(daily_freq_bands,1); freq_idx_day = 1:daily_num_bands;
            case {'multi_processor'}
                multi_proc = varargin{in_idx+1};
                if multi_proc == true
                    num_proc = varargin{in_idx+2};
                else
                	num_proc = 0;
                end
            case {'batch_processing'}
                batch_processing = true;
            case {'testing'}
                filename_add = '_TEST';
            case {'save_log'}
                save_log = true;
        end
    end
end

%% Initialize Multi-Processor Computing
if multi_proc
    if isempty(gcp('nocreate'))
        % create a local cluster object
        pc = parcluster('local');
        
        %explicitly set the JobStorageLocation to the temp directory that was
        %created in sbatch script (for batch processing only)
        if batch_processing
            pc.JobStorageLocation = pc_jobstorloc;
        end
        
        % start the matlabpool with num_proc workers
        parpool(pc, num_proc);
    end
end

%Pre-allocate log array
if save_log
    complete_log = cell(length(VarIndices)*length(modelArray),num_exps*length(seasonArray));
end

%% Calculate Variability
for i=1:length(VarIndices)
    %Define file variable identifiers
    [~,filevarFN,freq,~,~,~,~,~] = var_chars(VarIndices(i));
    if regexp(freq,'.[Dd]ay.'); freq_def = 1; year_length = 365; else freq_def = 2; year_length = 12; end %(To get correct frequency - day/month - info from arrays above)
    
    for j=1:num_models
        %Process Data by Experiment
        for experiment = 1:num_exps
            
            %Get file names/characteristics of files of inputted var/exp/model
            [~,~,~,~,num_conventions] = name_chars(modelArray{j},expArray(experiment),filevarFN,freq,'directory',[various_defaults.proc_data_dir,modelArray{j},various_defaults.season_dir]);
            
            if num_conventions == 0
                warningMsg = ['Files for ',modelArray{j},' ',expArray{experiment},' all seasons ',filevarFN,freq,' were not processed! (missing file)'];
                warning(warningMsg);
                if save_log
                    save_log_wrn = ['****',filevarFN,freq,modelArray{j},' ',expArray{experiment},' for all seasons Not Complete**** (missing file)'];
                    [complete_log{(i-1)*length(modelArray)+j,:}] = deal(save_log_wrn);
                end
            end
            
            %% Process Data
            for conv = 1:num_conventions
                [~,run,strtyr,endyr,~] = name_chars(modelArray{j},expArray(experiment),filevarFN,freq,'use_convention',conv,'directory',[various_defaults.proc_data_dir,modelArray{j},various_defaults.season_dir]);
                
                for season = 1:length(seasonArray)
                    try
                        initial_vars = who;
                        season_length = seasonLength{freq_def,season};
                        
                        %Determine _DATA file length and set up season
                        %calcs accordingly (depending on the years desired
                        %for each season in prefTimeframe)
                        if endyr{1}-strtyr{1}+1==prefTimeframe{season} %if _DATA file has preferred length, 'standard' length determined by convention for each season
                            %if season crosses over a year boundary, total
                            %number of seasons ("years") will be one less
                            %than years in file
                            if (regexp(freq,'.[Dd]ay.') && seasonStart{freq_def,season}+season_length>365) || (~regexp(freq,'.[Dd]ay.') && seasonStart{freq_def,season}+season_length>12)
                                years = prefTimeframe{season}-1; %years = how many seasons in time series
                            else %if season doesn't go over a year boundary
                                years = prefTimeframe{season};
                            end
                            %years = 30;
                            startbias = 0;
                        elseif endyr{1}-strtyr{1}+1>prefTimeframe{season} %if _DATA file is longer than needed for specific season (taking last [xx] years wanted)
                            if (regexp(freq,'.[Dd]ay.') && seasonStart{freq_def,season}+season_length>365) || (~regexp(freq,'.[Dd]ay.') && seasonStart{freq_def,season}+season_length>12)
                                years = prefTimeframe{season}-1; %years = how many seasons in time series
                            else %if season doesn't go over a year boundary
                                years = prefTimeframe{season};
                            end
                            startbias = endyr{1}-strtyr{1}+1-prefTimeframe{season};
                        else %if _DATA file is shorter than preferred
                            %if season crosses a year boundary
                            if (regexp(freq,'.[Dd]ay.') && seasonStart{freq_def,season}+season_length>365) || (~regexp(freq,'.[Dd]ay.') && seasonStart{freq_def,season}+season_length>12)
                                years = endyr{1}-strtyr{1};
                            else %if season doesn't go over a year boundary
                                years = endyr{1}-strtyr{1}+1;
                            end
                            startbias = 0;
                        end
                        strtyr_seas = strtyr{1}+startbias;
                        
                        %Define final filenames for existence testing and
                        %saving
                        filename1SD = [various_defaults.proc_data_dir,modelArray{j},folder,filevarFN,freq,modelArray{j},'_',expArray{experiment},'_',run{1},num2str(strtyr_seas),'-',num2str(endyr{1}),'_sqrtPower_',seasonArray{season}];
                        filename1Mean = [various_defaults.proc_data_dir,modelArray{j},folder,filevarFN,freq,modelArray{j},'_',expArray{experiment},'_',run{1},num2str(strtyr_seas),'-',num2str(endyr{1}),'_LocalMeans_',seasonArray{season}];
                        
                        %Decide whether calculations or replacements are needed
                        if replace_files %If replace_files, calculate and fully replace file
                            calc = true; append_calcs = false; calc_xf = true;
                        else %If ~replace_files, find out whether the file exists, and if so, whether the wanted bands exist
                            if exist([filename1SD,'.mat'],'file')
                                load(filename1SD);
                                if regexp(freq,'.[Mm]on.') %Find out which frequency bands already have completed calculations
                                    bands_exist = cell(monthly_num_bands,1);
                                    for band_test = 1:monthly_num_bands
                                        bands_exist{band_test} = structfind(StdDevs,'Size',monthly_freq_bands(band_test,:));
                                    end
                                else
                                    bands_exist = cell(daily_num_bands,1);
                                    for band_test = 1:daily_num_bands
                                        bands_exist{band_test} = structfind(StdDevs,'Size',daily_freq_bands(band_test,:));
                                    end
                                end
                                
                                %Find if interannual frequencies already calculated (procedure slightly different, so not in the same num_bands loop as the rest)
                                if structfind(StdDevs,'ID','XF') || structfind(StdDevs,'ID','XFavgDev'); calc_xf = false; 
                                else calc_xf = true; end
                                
                                if any(cellfun('isempty',bands_exist)) %If not all wanted bands exist, calculated, but add to file instead of replacing it
                                    calc = true; append_calcs = true;
                                    if ~all(cellfun('isempty',bands_exist)) %Subset to only needed bands
                                        if ~regexp(freq,'.[Dd]ay.');[~,idx_tmp] = ismember(monthly_freq_bands,cell2mat({StdDevs.Size}'),'rows'); freq_idx_month = find(~idx_tmp); clear idx_tmp
                                        else [~,idx_tmp] = ismember(daily_freq_bands,cell2mat({StdDevs.Size}'),'rows'); freq_idx_day = find(~idx_tmp); clear idx_tmp
                                        end
                                    end
                                else %If wanted bands all exist, no need to do the thing
                                    calc = false;
                                end
                            else %If file doesn't exist, calculate and create file
                                calc = true; append_calcs = false; calc_xf = true;
                            end
                        end
                        
                        if calc
                            filename = [various_defaults.proc_data_dir,modelArray{j},various_defaults.season_dir,filevarFN,freq,modelArray{j},'_',expArray{experiment},'_',run{1},num2str(strtyr{1}),'-',num2str(endyr{1}),'_DATA'];
                            load(filename);
                            
                            %Get length of dimensions
                            nlon=length(lon);
                            nlat=length(lat);
                            rge=size(Raw,3);
                            
                            %Pre-allocate FFT'ed and Means Matrices (for computing optimization)
                            AllPowers = zeros(nlon,nlat,years,season_length);
                            AllPowersMean = zeros(nlon,nlat,season_length);
                            LocalMeans_season = zeros(nlon,nlat,years); %(detrended/deseasonalized, for interannual var calcs)
                            
                            %Get index of first day of each season
                            seas_init = (((1:years)'-1+startbias)*year_length+seasonStart{freq_def,season});
                            seas_idxs = seas_init; %seas_idxs = years x length_season array of indices of all months needed
                            %Add indices of rest of each season (weird
                            %solution below because can't add (1:a) to each
                            %member of array otherwise)
                            for t = 1:season_length-1
                                seas_idxs = [seas_idxs seas_init+t]; %#ok<AGROW>
                            end
                            clear seas_init
                            %Get seasonal means of Raw data
                            LocalMean = mean(Raw(:,:,seas_idxs),3); %#ok<NASGU> %Local means over all seasons
                            
                            if regexp(freq,'.[Dd]ay.') %(data is daily)
                                freq_idx = freq_idx_day;
                                freq_bands = daily_freq_bands(freq_idx,:); num_bands = size(freq_bands,1);
                                freq_band_names = daily_freq_band_names(freq_idx);
                                freq_band_ID = daily_freq_band_ID(freq_idx);
                                XF_size = [365 Inf];
                                
                                %Calculate matrix for deseasonalization
                                X0 = ones(rge,1);
                                X_c = cos((2*pi/365)*(0:rge-1)'*(1:num_harmonics));
                                X_s = sin((2*pi/365)*(0:rge-1)'*(1:num_harmonics));
                                X = [X0 X_c X_s];
                                
                                %% Detrend and Deseasonalize
                                parfor (l=1:nlon,num_proc)
                                    %Detrend
                                    Y = detrend(squeeze(Raw(l,:,:))')';
                                    %Deseasonalize
                                    betahat = (X'*X)\(X'*Y(:,:)');
                                    Raw(l,:,:) = Y - betahat'*X';
                                end
                                clear Y
                            else %(data is monthly)
                                freq_idx = freq_idx_month;
                                freq_bands = monthly_freq_bands(freq_idx,:); num_bands = size(freq_bands,1);
                                freq_band_names = monthly_freq_band_names(freq_idx);
                                freq_band_ID = monthly_freq_band_ID(freq_idx);
                                XF_size = [12 Inf];
                                %Pre-allocate FFT'ed and Means Matrices (for computing optimization)
                                MonAvgs = zeros(nlon,nlat,12);
                                
                                %% Detrend and deseasonalize
                                for l=1:nlon
                                    %Detrend
                                    Y = detrend(squeeze(Raw(l,:,:))')';
                                    %Take Average of each Month over the Time Series
                                    for n=1:12
                                        MonAvgs(l,:,n) = mean(Y(:,n:12:rge),2);
                                    end
                                    %Remove Average of each Month from Time Series to Deseasonalize
                                    a = repmat(squeeze(MonAvgs(l,:,:)),1,rge/12);
                                    Raw(l,:,:) = Y - a;
                                end
                                clear Y
                            end
                            %% Split up by season
                            Raw_seasons = zeros(nlon,nlat,years,season_length);
                            for t = 1:years
                                Raw_seasons(:,:,t,:) = Raw(:,:,seas_idxs(t,:));
                            end
                            clear Raw
                            %% Get Power Spectrum and Mean
                            parfor (l=1:nlon, num_proc)
                                for t=1:years
                                    %Mean of each (deseasonalized) season
                                    LocalMeans_season(l,:,t) = mean(squeeze(Raw_seasons(l,:,t,:)),2);
                                    %DFT each season separately
                                    Y_c2 = season_length.*abs(fft(squeeze(Raw_seasons(l,:,t,:)),season_length,2)./season_length).^2;
                                    AllPowers(l,:,t,:) = Y_c2;
                                end
                            end
                            clear Raw_seasons Y_c2
                            %Get average of all seasons
                            for l = 1:nlon
                                AllPowersMean(l,:,:) = mean(AllPowers(l,:,:,:),3);
                            end
                            clear AllPowers
                            
                            %Put frequencies in order
                            nfreqs = size(AllPowersMean,3);
                            Ks = (0:(nfreqs/2))./nfreqs;
                            period = 1./Ks; %period
                            
                            %Calculate and store band-separated standard deviations
                            freq_bands_idxs = zeros(size(freq_bands));
                            for band = 1:num_bands
                                %Get indices in frequency space of desired
                                %frequency bands
                                [~,freq_bands_idxs(band,1)] = min(abs(freq_bands(band,1)-period));
                                [~,freq_bands_idxs(band,2)] = min(abs(freq_bands(band,2)-period));
                                %If sequential, make sure a frequency is not
                                %double-counted
                                if band>1 && freq_bands_idxs(band,1) == freq_bands_idxs(band-1,2)
                                    freq_bands_idxs(band,1) = freq_bands_idxs(band,1)-1;
                                end
                                %Preallocate matrix
                                band_stddev = zeros(nlon,nlat);
                                %Sum densities, take square roots to make
                                %standard deviation from variance
                                for l = 1:nlon
                                    for m = 1:nlat
                                        band_stddev(l,m) = ((sum(AllPowersMean(l,m,freq_bands_idxs(band,2):freq_bands_idxs(band,1))))/(nfreqs/2))^0.5;
                                    end
                                end
                                %Save local band-separated standard deviation with corresponding band information
                                StdDevs_tmp(band).Name = freq_band_names{band}; %#ok<AGROW>
                                StdDevs_tmp(band).ID = freq_band_ID{band}; %#ok<AGROW>
                                StdDevs_tmp(band).Size = freq_bands(band,:); %#ok<AGROW>
                                StdDevs_tmp(band).Data = band_stddev; %#ok<AGROW>
                            end
                            
                            %Calculate and store interannual standard
                            %deviations (the interannual seasonal
                            %variability), in the same struct
                            if calc_xf
                                StdDevs_tmp(num_bands+1).Data = std(LocalMeans_season,0,3);
                                StdDevs_tmp(num_bands+1).ID = 'XF';
                                StdDevs_tmp(num_bands+1).Size = XF_size; %by convention
                                StdDevs_tmp(num_bands+1).Name = '> 1 year';
                            end
                            
                            ExpMsg = ['Calculations for ',modelArray{j},' ',expArray{experiment},' ',seasonArray{season},' ',filevarFN,freq,' Convention ',num2str(conv),' Complete'];
                            disp(ExpMsg)
                            
                            %Add calculated bands to existing ones, if desired
                            if append_calcs
                                StdDevs = [StdDevs';StdDevs_tmp']'; clear StdDevs_tmp;
                            else
                                StdDevs = StdDevs_tmp; clear StdDevs_tmp
                            end
                            
                            %% Save Data
                            filename1SD = [filename1SD,filename_add,'.mat']; %#ok<AGROW>
                            filename1Mean = [filename1Mean,filename_add,'.mat']; %#ok<AGROW>
                            save(filename1SD,'-v7.3','StdDevs','lat','lon');
                            save(filename1Mean,'-v7.3','LocalMean','lat','lon');
                            clear lat lon
                            
                            disp([filename1SD,' saved']);
                            disp([filename1Mean,' saved']);
                            %Store success message in log
                            if save_log
                                save_log_msg = [filevarFN,freq,modelArray{j},' ',expArray{experiment},' ',seasonArray{season},' ',num2str(conv),' Complete'];
                                complete_log{(i-1)*length(modelArray)+j,(experiment-1)*length(seasonArray)+season} = save_log_msg;
                            end
                        else
                            disp(['Files for ',modelArray{j},' ',expArray{experiment},' ',filevarFN,freq,num2str(strtyr{1}),'-',num2str(endyr{1}),' already exist, were not replaced']);
                            if save_log
                                save_log_msg = [filevarFN,freq,modelArray{j},' ',expArray{experiment},' ',seasonArray{season},' ',num2str(conv),' already exists'];
                                complete_log{(i-1)*length(modelArray)+j,(experiment-1)*length(seasonArray)+season} = save_log_msg;
                            end
                        end
                        endMsg = ['Files for ',modelArray{j},' ',expArray{experiment},' ',seasonArray{season},' ',num2str(conv),' ',filevarFN,freq,' Processed!'];
                        disp(endMsg);
                        %Clear variables not needed between experiments
                        clearvars('-except',initial_vars{:});
                    catch ME %catch for issues with specific var-exp-model-seas
                        disp(ME);
                        try %#ok<TRYNC>
                            disp(ME.message);
                            disp(ME.stack);
                            try %#ok<TRYNC>
                                disp(ME.stack.line);
                            end
                        end
                        warningMsg = ['Files for ',modelArray{j},' ',expArray{experiment},' ',seasonArray{season},' ',num2str(conv),' ',filevarFN,freq,' were not processed!'];
                        warning(warningMsg);
                        %Store warning in log
                        if save_log
                            save_log_wrn = ['****',filevarFN,freq,modelArray{j},' ',expArray{experiment},' ',seasonArray{season},' ',num2str(conv),' Not Complete****'];
                            complete_log{(i-1)*length(modelArray)+j,(experiment-1)*length(seasonArray)+season} = save_log_wrn;
                        end
                        %Clear variables not needed between experiments
                        clearvars('-except',initial_vars{:});
                    end
                end
            end
        end
    end
    disp(['All Files for ',filevarFN,freq,' Complete!']);
end
%Export log as table
if save_log
    complete_log = complete_log';
    complete_log = complete_log(:);
    export_log = cell2table(complete_log,'VariableNames',{'VariabilitySeasonCalcs_attempted_executions'});
    writetable(export_log,[various_defaults.log_dir,'Variability_seasons_log_',num2str(startTimestamp(1)),'-',num2str(startTimestamp(2)),'-',num2str(startTimestamp(3)),'-',num2str(startTimestamp(4)),'-',num2str(startTimestamp(5)),'-',num2str(startTimestamp(6)),'.txt']);
end
toc
end

