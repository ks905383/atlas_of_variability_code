function Variability(VarIndices,modelArray,varargin)
% VARIABILITY   calculate variability by frequency band.
%
%   VARIABILITY(VarIndices,modelArray) calculates the contribution of
%   signals in various frequency bands to the standard deviation of
%   variables from models given by [VarIndices] and [modelArray]. All data
%   is first deseasonalized - by a harmonic filter for daily data and by
%   removing the average of all of each month from that month for monthly
%   data. Standard deviations are calculated by taking the square root of
%   the sums of the power spectral density scaled by the length of the time
%   series up to set frequencies. The local mean of each variable is also
%   saved. The program outputs two files for independent [_DATA.mat] files
%   processed - one file [_sqrtPower] with a structure array StdDev
%   containing the data separated by frequency band, word description of
%   each frequency band, and the frequency band sizes/endpoints; in
%   addition to (outside of the struct) the lat and lon variables from the
%   original file - and one file [_LocalMeans] with the local mean of the
%   variable in addition to the lat and lon variables. Code catches
%   failures with warnings at the individual variable - model - experiment
%   level.
%
%   Default frequency bands (by default, 30 years or largest possible
%   subset taken from data):
%       ID      Daily       Monthly     Notes
%       HF      1-5         -           saved as [0,5] by convention
%       MF      5-30        -
%       LF      30-365      1-12
%       XF      > 1 yr      > 1 yr      saved as [12,Inf] or [365,Inf]
%       Full    all         all         saved as [0,Inf] by convention
%		2to15	2-15		- 			
%		15to90	15-90		- 
%		XXF		- 			12-120
%		XXXF	-			> 10 yrs	saved as [120,Inf] by convention
%
%
%   Sample Output for 2-band + overall monthly [_sqrtPower] file:
%   StdDevs.Name = {'< 1 year';'> 1 year';'all bands'};
%   StdDevs.ID = {'LF';'XF';'Full'};
%   StdDevs.Size = {0,12;12,Inf;0,Inf};
%   StdDevs.Data = {[nlon x nlat array];[nlon x nlat array;,...
%                       [nlon x nlat array]};
%   lat = [nlat x 1 array]
%   lon = [nlon x 1 array]
%
%   Data requirements: [_DATA.mat] files for each variable, experiment,
%                      model
%   Function requirements (on path): name_chars, var_chars
%
%   VARIABILITY(...,'[flag]',[params],...) modify program run as below:
%       'experiments',[cell]- manually set experiments to process (default
%                             {'rcp85','historical','piControl'})
%       'num_harmonics',[#] - change the number of harmonics removed in
%                             deseasonalization process for daily data
%       'remove_locavg',[log]   - set whether to smooth fourier transform
%                                 at frequency endpoints (def: true)
%       'num_splits',[#]    - set the number of subsets of data to be
%                             loaded at a time (def: 1 --> entire _DATA
%                             file is loaded). If dealing with _DATA files
%                             of significant size (> 55K pixels, > 30
%                             years, etc.), setting a high num_splits can
%                             drastically reduce memory use at the expense
%                             of computing speed (calls to matfile are
%                             pretty slow). Maximum memory used in
%                             workspace: almost the entirety of memory use
%                             is the ceil(nlon/num_split) x nlat x time
%                             array of Raw data, followed by time x 25
%                             deseasonalization (if daily) array and nlon x
%                             nlat pre-allocated arrays of output data.
%       'freq_bands',[names],[lims],[ID]- set new endpoints for frequency
%                                         bands. [names] is a num_freq x 1
%                                         cell array of strings giving a
%                                         short description of the
%                                         frequency band. [lims] is a
%                                         num_freq x 2 numerical array
%                                         giving the start (first column)
%                                         and end (second column) times (in
%                                         period) of the frequency band.
%                                         [ID] is a num_freq x 1 cell array
%                                         of strings giving a short, easily
%                                         searchable ID for each frequency
%                                         band. (default for daily data
%                                         i.e. is [names] = {'< 5 days','5
%                                         - 30 days','30 - 365 days','> 1
%                                         year','all bands'}, [lims] =
%                                         [0,5;5,30;30,365;365,Inf;0,Inf],
%                                         and [ID] =
%                                         {'HF','MF','LF','XF','Full'});
%       'multi_processor',[log],([#])   - set whether to use multiple
%                                         processors and set number of [#]
%                                         processors to be used if true 
%										  (def: true,12)
%       'batch_processing',[log]        - sets job storage location
%                                         explicitly at
%                                         /tmp/kschwarzwald/[SLURMJOBID]
%                                         (def: false; only useful in batch
%                                         processing)
%       'replace',[log]     - set whether to replace existing files (def:
%                             true)
%       'save_log'          - exports log of attempted saves with
%                             success/failure status under
%                             [log_dir]/Variability_log_[starttime].txt
%       'testing'           - adds '_TEST' to end of filename
%
%   Saving convention:
%   [filevar]_[freq]_[model]_[exp]_[run]_[strtyr]_[endyr]_sqrtPower.mat
%   [filevar]_[freq]_[model]_[exp]_[run]_[strtyr]_[endyr]_LocalMeans.mat
%   NOTE: VARIABILITY() will create a model subdirectory if it does not yet
%   exist. 
%
%   RUNTIME:
%       - 2.5 seconds for a 29x96x10950 chunk (29 x 96 ffts over 30
%         years of daily data, vectorized by longitude) + auxiliary time at
%         start and end to setup, save; timed using 12 cpus on an RCC
%         bigmem2 node.
%       - ~75 seconds (longer for first few due to startup costs) to load a
%         (1-8) x 48 x 36500 array subset from a .mat file (i.e., using
%         1000-year daily CCSM3 data, split up for individual longitude
%         bands) (load from matfile is somewhat independent of how much is
%         being loaded due to overhead).
%       - NB: parallel computing functionality will likely not bring
%         maximized time efficiency gains if the number of longitude bands
%         per split is fewer than the number of processors used (the base
%         array sent to each parallel worker (after splitting by longitude
%         bands using [num_splits]) is [1 longitude band] x [num latitudes]
%         x [num timesteps] for vectorization purposes).
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also VARIABILITY_CI, VARIABILITY_SEASONS
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 11/26/2017

%% 1 Misc Setup
% Allow for string input (as opposed to cell) for model name
if ~isa(modelArray,'cell'); modelArray = {modelArray}; end

% Get timestamp (to identify logs)
format shortG
startTimestamp = clock;
start_msg = ['\nVariability calculations run, \nstarted at ',...
    num2str(startTimestamp(4)),':',num2str(startTimestamp(5)),':',num2str(fix(startTimestamp(6))),...
    ' on ',num2str(startTimestamp(2)),'/',num2str(startTimestamp(3)),'/',num2str(startTimestamp(1)),...
    ',\nfor variable indices [',num2str(VarIndices),'] \nand models ',sprintf('%s ',modelArray{:}),'\nGood Luck! :) \n ',quoteme,' \n \n \n'];
fprintf(start_msg)

%% 2 Set default values and set up program
%Get stored defaults
various_defaults = matfile('various_defaults.mat');
%Set other default options
replace_files = true;
num_models = length(modelArray);
num_harmonics=12; remove_locavg = true;
expArray = {'rcp85','historical','piControl'};
num_exps = length(expArray);
multi_proc = true; num_proc = 12; batch_processing = false;
%pc_jobstorloc = strcat('/tmp/kschwarzwald/', getenv('SLURM_JOB_ID'));
pc_jobstorloc = strcat(various_defaults.tmp_dir,randi(10000));
save_log = false;
filename_add = [];
num_splits = 1;

file_end_save = [];
file_end = '*DATA.mat';

%Frequency band determinants default settings (taken from various_defaults)
cust_freqs = false;
freq_band_setup = various_defaults.freq_band_setup;

%% 3 Set behavior of optional function flags
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
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
            case {'remove_locavg'}
                remove_locavg = varargin{in_idx+1};
            case {'freq_bands'}
                freq_band_names = varargin{in_idx+1}; varargin{in_idx+1} = 0;
                freq_bands = varargin{in_idx+2}; varargin{in_idx+2} = 0; if ~isa(freq_bands,'numeric'); error('VARIABILITY:NotArray','freq band sizes must be a numerical array'); end
                freq_band_ID = varargin{in_idx+3}; varargin{in_idx+3} = 0;
                num_bands = size(freq_bands,1); freq_idx = 1:num_bands;
                cust_freqs = true;
            case {'num_splits'}
                num_splits = varargin{in_idx+1};
            case {'multi_processor'}
                multi_proc = varargin{in_idx+1};
                if multi_proc == true
                    num_proc = varargin{in_idx+2};
                else
                	num_proc = 0;
                end
            case {'batch_processing'}
                batch_processing = true;
            case {'save_log'}
                save_log = true;
            case {'testing'}
                filename_add = '_TEST';
            case {'file_end'}
                file_end = varargin{in_idx+1};
                file_end_save = varargin{in_idx+1};
        end
    end
end

%% 4 Initialize Multi-Processor Computing
if multi_proc
    if isempty(gcp('nocreate'))
        % create a local cluster object
        pc = parcluster('local');
        
        %explicitly set the JobStorageLocation to the temp directory that was
        %created in sbatch script (for batch processing only, this is to
        %ensure that multiple parallel initializations don't overlap and
        %modify each others' files - see
        %https://rcc.uchicago.edu/docs/software/environments/matlab/ for
        %more info). 
        if batch_processing
            pc.JobStorageLocation = pc_jobstorloc;
        end
        
        % start the matlabpool with num_proc workers
        parpool(pc, num_proc);
    end
end

%Pre-allocate log array
if save_log
    complete_log = cell(1,length(VarIndices)*length(modelArray));
end

%% 5 Calculate Variability
for var_idx=1:length(VarIndices)
    %Define file variable identifiers
    [~,filevarFN,freq] = var_chars(VarIndices(var_idx));
    
    %Set up freq bands and deseasonalization, depending on frequency (if
    %not custom)
    if ~cust_freqs
        if regexp(freq,'.[Dd]ay.')
            freq_char_id = find(~cellfun(@isempty,strfind(freq_band_setup(:,1),'day')));
        elseif regexp(freq,'.[Mm]on.')
            freq_char_id = find(~cellfun(@isempty,strfind(freq_band_setup(:,1),'month')));
        elseif regexp(freq,'.yrr.')
            freq_char_id = find(~cellfun(@isempty,strfind(freq_band_setup(:,1),'year')));
        else
            freq_char_id = find(~cellfun(@isempty,strfind(freq_band_setup(:,1),'3hour')));
            warning('Variability:NotSuppFreq',['Processing for the sampling frequency ',freq,...
                ' is not yet validated - only daily, ',...
                'monthly, or yearly data are currently fully supported.'])
        end
        
        freq_bands = freq_band_setup{freq_char_id,3};
        freq_band_names = freq_band_setup{freq_char_id,2};
        freq_band_ID = freq_band_setup{freq_char_id,4};
        num_bands = size(freq_bands,1);
    end
    
    for model_idx=1:num_models
        model = modelArray{model_idx};
        %Process Data by Experiment
        for exp_id = 1:num_exps
            initial_vars = who;
            experiment = expArray{exp_id};
            try
                %% 5.1 Get file names/characteristics of files of inputted var/exp/model
                [~,~,~,~,num_conventions] = name_chars(model,expArray(exp_id),filevarFN,freq,'file_end',file_end);
                if num_conventions == 0
                    error('VARIABILITY:MissingData',...
                        ['name_chars found no files for the combination, model: ',...
                        model,', experiment: ',expArray{exp_id},...
                        ', variable: ',filevarFN,', and frequency: ',freq])
                end
                %% 5.2 Process Data
                for conv = 1:num_conventions
                    %% 5.2.1 Define filenames
                    [~,run,strtyr,endyr,~] = name_chars(model,expArray(exp_id),filevarFN,freq,'use_convention',conv,'file_end',file_end);
                    
                    %Define final filenames (any filename_add such as '_TEST' is added at the end to not interfere with existence testing)
                    filename1SD = [various_defaults.proc_data_dir,model,'/',filevarFN,freq,model,'_',experiment,'_',run{1},num2str(strtyr{1}),'-',num2str(endyr{1}),'_sqrtPower',file_end_save];
                    filename1Mean = [various_defaults.proc_data_dir,model,'/',filevarFN,freq,model,'_',experiment,'_',run{1},num2str(strtyr{1}),'-',num2str(endyr{1}),'_LocalMeans',file_end_save];
                    %Make save directory, if it does not yet exist
                    if ~exist([various_defaults.proc_data_dir,model,'/'],'file')
                        mkdir([various_defaults.proc_data_dir,model,'/'],'file')
                    end
                    
                    %% 5.2.2 Decide whether calculations or replacements are needed
                    if replace_files %If replace_files, calculate and fully replace file
                        calc = true; replace = true;
                    else %If ~replace_files, find out whether the file exists, and if so, whether the wanted bands exist
                        if exist([filename1SD,'.mat'],'file')
                            load(filename1SD);
                            bands_exist = cell(num_bands,1);
                            for band_test = 1:num_bands
                                bands_exist{band_test} = structfind(StdDevs,'Size',freq_bands(band_test,:));
                            end
                            if any(cellfun('isempty',bands_exist)) %If not all wanted bands exist, calculated, but add to file instead of replacing it
                                calc = true; replace = false;
                                if ~all(cellfun('isempty',bands_exist)) %Subset to only needed bands
                                    [~,idx_tmp] = ismember(freq_bands,cell2mat({StdDevs.Size}'),'rows'); freq_idx = find(~idx_tmp); clear idx_tmp
                                end
                                if freq_idx<=size(freq_bands,1)
                                    freq_bands = freq_bands(freq_idx,:); num_bands = size(freq_bands,1);
                                    freq_band_names = freq_band_names(freq_idx);
                                    freq_band_ID = freq_band_ID(freq_idx);
                                else
                                    error('VARIABILITY:MultConvBug',...
                                        ['This is a bug caused by replace_files = false and num_convs>1,',...
                                        ' which causes freq_bands to potentially change size between the conv for loop runs (causing a "index exceeds matrix dimensions"). ',...
                                        'In the short run, rerunning the program will fix the issue, but a longer-term fix should be attempted (likely by making a _tmp version of all freq_band variables to draw on from inside the for loop)'])
                                end
                            else %If wanted bands all exist, no need to do the thing
                                calc = false;
                            end
                        else %If file doesn't exist, calculate and create file
                            calc = true; replace = true;
                        end
                    end
                    
                    %% 5.2.3 Calculate variability, if needed
                    if calc
                        %% 5.2.3.1 Load data and pre-allocate arrays
                        filename = [various_defaults.raw_data_dir,model,'/',filevarFN,freq,model,'_',experiment,'_',run{1},num2str(strtyr{1}),'-',num2str(endyr{1}),'_DATA'];
                        %Get matfile object of _DATA
                        load_file = matfile([filename,'.mat']);
                        %Get dimensions of data
                        load_file_vars = whos(load_file);
                        if isempty(structfind(load_file_vars,'name','Raw'))
                            error('VARIABILITY:NoRawVar',['No variable of the name "Raw" found in the file ',filename,'.mat. I would suggest running /Code/File_Utilities/filevar_Raw_fixes.m; the likely cause is an old file system convention.']) 
                        end
                        nlon = load_file_vars(structfind(load_file_vars,'name','Raw')).size(1);
                        nlat = load_file_vars(structfind(load_file_vars,'name','Raw')).size(2);
                        rge = load_file_vars(structfind(load_file_vars,'name','Raw')).size(3);
                        
                        %Pre-allocate FFT'ed and Means Matrices (for computing optimization)
                        LocalMean = zeros(nlon,nlat); StdDev = zeros(nlon,nlat);
                        
                        if regexp(freq,'.[Dd]ay.')
                            %Calculate Matrix for Deseasonalization
                            X0 = ones(rge,1);
                            X_c = cos((2*pi/365)*(0:rge-1)'*(1:num_harmonics));
                            X_s = sin((2*pi/365)*(0:rge-1)'*(1:num_harmonics));
                            X = [X0 X_c X_s]; clear X0 X_c X_s
                        elseif regexp(freq,'.3hr.')
                            %Calculate Matrix for Deseasonalization
                            X0 = ones(rge,1);
                            X_c = cos((2*pi/(365*8))*(0:rge-1)'*(1:num_harmonics));
                            X_s = sin((2*pi/(365*8))*(0:rge-1)'*(1:num_harmonics));
                            X = [X0 X_c X_s]; clear X0 X_c X_s
                        elseif regexp(freq,'.[Mm]on.')
                            %Pre-allocate Means matrix
                            MonAvgs = zeros(nlon,nlat,12);
                        end
                        
                        %Put frequencies in order
                        Ks = (0:(rge/2))./rge;
                        period = 1./Ks; %Calculate periods
                        
                        %Set indices of frequency bands to calculate
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
                            band_stddev = zeros(nlon,nlat,num_bands);
                        end
                        
                        %% 5.2.3.2 Process variable by longitude split
                        for data_splits = 1:num_splits
                            %% 5.2.3.2a Set start and end longitude values
                            if num_splits == 1
                                lon_i = 1;
                                lon_f = nlon;
                            else
                                lon_i = (data_splits-1)*ceil(nlon/num_splits)+1;
                                lon_f = (data_splits)*ceil(nlon/num_splits);
                                if lon_f > nlon; lon_f = nlon;
                                end
                            end
                            
                            disp(['Begin loading lon: ',num2str(lon_i),' - ',num2str(lon_f)])
                            tic
                            %Load variable
                            Raw = load_file.Raw(lon_i:lon_f,:,:);
                            %Remove crazy outliers 
                            Raw(Raw>10e30) = 0; 
                            %Pre-allocate tmp means
                            LocalMean_tmp = zeros(length(lon_i:lon_f),nlat);
                            toc
                            
                            disp(['Begin calculations for lon: ',num2str(lon_i),' - ',num2str(lon_f)])
                            tic
                            
                            %% 5.2.3.2b Detrend and Deseasonalize, depending on frequency
                            if ~isempty(regexp(freq,'.[Dd]ay.', 'once')) || ~isempty(regexp(freq,'.3hr.','once')) %(data is daily)
                                %% Detrend and Deseasonalize
                                parfor (lon_idx=1:length(lon_i:lon_f),num_proc)
                                    %Calculate Local Means (tmp file
                                    %because of parfor requirements)
                                    LocalMean_tmp(lon_idx,:) = mean(Raw(lon_idx,:,:),3);
                                    %Detrend
                                    Y = detrend(squeeze(Raw(lon_idx,:,:))')'; %MUST ' INSIDE DETREND TO ENSURE DETREND ACROSS TIME, NOT LATS
                                    %Y = detrend(squeeze(Raw(lon_idx,:,:))');
                                    %Deseasonalize
                                    betahat = (X'*X)\(X'*Y');
                                    %Replace original timeseries with
                                    %detrended, deseasonalized timeseries
                                    Raw(lon_idx,:,:) = Y - betahat'*X';
                                end
                                %Add tmp means to full means file
                                LocalMean(lon_i:lon_f,:) = LocalMean_tmp;
                                clear Y LocalMean_tmp
                            elseif ~isempty(regexp(freq,'.[Mm]on.','once')) %(data is monthly)
                                %% Detrend and deseasonalize
                                for lon_idx=1:length(lon_i:lon_f)
                                    %Calculate Local Means
                                    LocalMean(lon_idx+lon_i-1,:) = mean(Raw(lon_idx,:,:),3);
                                    %Detrend
                                    Y = detrend(squeeze(Raw(lon_idx,:,:))')';
                                    %Take Average of each Month over the Time Series
                                    for n=1:12
                                        MonAvgs(lon_idx,:,n) = mean(Y(:,n:12:rge),2);
                                    end
                                    %Remove Average of each Month from Time Series to Deseasonalize
                                    a = repmat(squeeze(MonAvgs(lon_idx,:,:)),1,rge/12);
                                    %Replace original timeseries with
                                    %detrended, deseasonalized timeseries
                                    Raw(lon_idx,:,:) = Y - a;
                                end
                                clear Y
                            elseif ~isempty(regexp(freq,'.yrr.','once')) %(data is yearly)
                                for lon_idx=1:length(lon_i:lon_f)
                                    %Calculate Local Means
                                    LocalMean(lon_idx+lon_i-1,:) = mean(Raw(lon_idx,:,:),3);
                                end
                            end
                            
                            %% 5.2.3.2c Calculate Standard Deviation
                            for lon_idx=1:length(lon_i:lon_f)
                                %Calculate Local Means
                                StdDev(lon_idx+lon_i-1,:) = std(Raw(lon_idx,:,:),0,3);
                            end
                            
                            %% 5.2.3.2d Calculate Power Spectrum
                            parfor (lon_idx = 1:length(lon_i:lon_f),num_proc)
                                %DFT
                                Y_c2 = rge.*abs(fft(squeeze(Raw(lon_idx,:,:)),rge,2)./rge).^2;
                                %Remove values at seasonal harmonics, replaces with local average (smoothing over gaps)
                                if ~isempty(regexp(freq,'.[Dd]ay.','once')) && remove_locavg
                                    harmonic_idxs = (1:num_harmonics)./365;
                                    for x=1:length(harmonic_idxs)
                                        [~,idx] = min(abs(harmonic_idxs(x)-Ks));
                                        Y_c2(:,idx) = (Y_c2(:,idx-1)+Y_c2(:,idx+1))/2;
                                    end
                                    Y_c2(:,1)=(Y_c2(:,2)+Y_c2(:,3))/2;
                                end
                                %Save power spectrum
                                Raw(lon_idx,:,:) = Y_c2;
                            end
                            clear Y_c2
                            
                            %% 5.2.3.2e Get variability from power spectrum
                            %Calculate and store band-separated standard deviations
                            for band = 1:num_bands
                                %Sum over spectral density to get variance,
                                %take square root to get standard deviation
                                for lon = lon_i:lon_f
                                    for lat = 1:nlat
                                        band_stddev(lon,lat,band) = ((sum(Raw(lon-lon_i+1,lat,freq_bands_idxs(band,2):freq_bands_idxs(band,1))))/(rge/2))^0.5;
                                    end
                                end
                                clear lon lat
                            end
                            toc
                            
                        end
                        
                        %% 5.2.3.3 Create and populate final data struct
                        %Save local band-separated standard deviations with
                        %corresponding band information
                        for band = 1:num_bands
                            StdDevs_tmp(band).Name = freq_band_names{band}; %#ok<AGROW>
                            StdDevs_tmp(band).ID = freq_band_ID{band}; %#ok<AGROW>
                            StdDevs_tmp(band).Size = freq_bands(band,:); %#ok<AGROW>
                            StdDevs_tmp(band).Data = squeeze(band_stddev(:,:,band)); %#ok<AGROW>
                        end
                        %Add means to same struct
                        StdDevs_tmp(num_bands+1).Name = 'mean';
                        StdDevs_tmp(num_bands+1).ID = 'mean';
                        StdDevs_tmp(num_bands+1).Size = [];
                        StdDevs_tmp(num_bands+1).Data = LocalMean;
                        %Add std devs to same struct
                        StdDevs_tmp(num_bands+2).Name = 'standard deviation';
                        StdDevs_tmp(num_bands+2).ID = 'std';
                        StdDevs_tmp(num_bands+2).Size = [];
                        StdDevs_tmp(num_bands+2).Data = StdDev;
                        ExpMsg = ['Calculations for ',model,' ',experiment,' ',filevarFN,freq,' Convention ',num2str(conv),' Complete'];
                        disp(ExpMsg)
                        
                        %% 5.2.3.4 Compile and save data
                        %Add calculated bands to existing ones, if desired
                        if ~replace
                            StdDevs = [StdDevs';StdDevs_tmp']'; clear StdDevs_tmp;
                        else
                            StdDevs = StdDevs_tmp; clear StdDevs_tmp
                        end
                        
                        %Get lat, lon data from _DATA file
                        lat = load_file.lat; %#ok<NASGU>
                        lon = load_file.lon; %#ok<NASGU>
                        
                        %Save Data
                        filename1SD = [filename1SD,filename_add,'.mat']; %#ok<AGROW>
                        filename1Mean = [filename1Mean,filename_add,'.mat']; %#ok<AGROW>
                        save(filename1SD,'-v7.3','StdDevs','lat','lon');
                        save(filename1Mean,'-v7.3','LocalMean','lat','lon');
                        
                        disp([filename1SD,' saved']);
                        disp([filename1Mean,' saved']);
                        %Store success message in log
                        if save_log
                            save_log_msg = [filevarFN,freq,model,' Complete'];
                            complete_log{(var_idx-1)*length(modelArray)+model_idx} = save_log_msg;
                        end
                        
                    else
                        disp(['Files for ',model,' ',experiment,' ',filevarFN,freq,num2str(strtyr{1}),'-',num2str(endyr{1}),' already exist, were not replaced']);
                        if save_log
                            save_log_msg = [filevarFN,freq,model,' already exists'];
                            complete_log{(var_idx-1)*length(modelArray)+model_idx} = save_log_msg;
                        end
                    end
                end
                
                %% 5.3 Display end message and clean up workspace
                endMsg = ['Files for ',model,' ',experiment,' ',filevarFN,freq,' Processed!'];
                disp(endMsg);
                %Clear variables not needed between experiments
                clearvars('-except',initial_vars{:});
                
            catch ME %catch for issues with specific experiment
                disp(ME);
                try %#ok<TRYNC>
                    disp(ME.message);
                    disp(ME.stack);
                     try %#ok<TRYNC>
                         for a = length(ME.stack)
                             disp(ME.stack(a).line);
                             disp(ME.stack(a).name);
                         end
                     end
                end
                warningMsg = ['Files for ',model,' ',experiment,' ',filevarFN,freq,' were not processed!'];
                warning(warningMsg);
                %Store warning in log
                if save_log
                    save_log_wrn = ['****',filevarFN,freq,model,' Not Complete****'];
                    complete_log{(var_idx-1)*length(modelArray)+model_idx} = save_log_wrn;
                end
                %Clear variables not needed between experiments
                clearvars('-except',initial_vars{:});
            end
        end
    end
    disp(['All Files for ',filevarFN,freq,' Complete!']);
end

%% 6 Export log as table, if desired
if save_log
    export_log = cell2table(complete_log','VariableNames',{'VariabilityCalcs_attempted_executions'});
    writetable(export_log,[various_defaults.log_dir,'/Variability_log_',num2str(startTimestamp(1)),'-',num2str(startTimestamp(2)),'-',num2str(startTimestamp(3)),'-',num2str(startTimestamp(4)),'-',num2str(startTimestamp(5)),'-',num2str(startTimestamp(6)),'.txt']);
end

end
