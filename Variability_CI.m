function Variability_CI(VarIndices,modelArray,varargin)
% VARIABILITY_CI   Find indices of not 'meaningful' variability ratios.
%
%   VARIABILITY_CI(VarIndices,modelArray) calculates the ratio of
%   variability between two experimental runs of a climate model by
%   frequency bands (and in total - stddev) and gives the index of pixels
%   not meaningfully different from 1. These pixels are defined as those
%   including 1 within a 2-sig confidence interval as determined using a
%   moving-block bootstrapping process. Variability calculations are
%   equivalent to those conducted in VARIABILITY. The program outputs
%   data in a file [_StdDevsCI.mat] with a structure array [StdDevsRatio]
%   containing data separated by frequency band, frequency band
%   sizes/endpoints, word description of each frequency band, lat and lon
%   indices of 'not meaningful' pixels, and the local standard deviation of
%   bootstrapped runs in addition to (outside of the struct) the lat and
%   lon variables from the original [_DATA] file used. In addition, an
%   included structure array [file_params] gives the start and end years
%   and run used for each experiment. 
%
%   If a [StdDevs] file does not yet exist (as calculated from
%   VARIABILITY), then one is created in the same format using VARIABILITY.
%
%   Default frequency bands (by default, 30 years or largest possible
%   subset taken from data):
%       ID      Daily       Monthly     Notes
%       HF      1-5         -           saved as [0,5] by convention
%       MF      5-30        - 
%       LF      30-365      1-12
%       XF      > 1 yr      > 1 yr      saved as [12,Inf] or [365,Inf]
%       2to15   < 15        -
%       15to90  15-90       -
%       Full    all         all         saved as [0,Inf] or [0 Inf] by
%                                       convention
%
%   Output of StdDevsRatio (for each frequency band [band]):
%   StdDevsRatio(band).Name = [char] (freq band description, i.e. '< 5
%                                     days')
%   StdDevsRatio(band).ID = [char] (Identifier, i.e. 'HF','LF',...)
%   StdDevsRatio(band).Size = [num] (freq band size, i.e. 5, 30, 365,...)
%   StdDevsRatio(band).Data = [nlon x nlat array] (Variability Ratios)
%   StdDevsRatio(band).lat_ciidx = [num_insig x 1 array] (lat indices)
%   StdDevsRatio(band).lon_ciidx = [num_insig x 1 array] (lon indices)
%   StdDevsRatio(band).std_runs = [nlon x nlat array] (std deviations of
%                                                  ratio of bootstrap runs) 
%
%   Data requirements: [_DATA.mat] files for each variable, experiment, and
%                      model, preferably _sqrtPower for all variables to
%                      process (call to Variability from inside this
%                      program is memory-intensive)
%   Function requirements (on path): name_chars, var_chars, and (from
%                                    community) structfind
%
%   Variability_CI(...,'[flag]',[params],...) modify program run as below:
%       'experiments',[cell]  - manually set experiments to process. Ratio
%                               is calculated [cell]{1} / [cell]{2}
%                               (default {'rcp85','piControl'})
%       'block_size',[#]      - set bootstrap block size [#] (def: 365 time
%                               units)
%       'nruns',[#]           - set number of bootstrap runs [#] (def: 1200)
%       'rand_seed',[#]       - set seed for random number generator [#]
%                               (def: 65987436)
%       'remove_locavg',[log] - set whether to replace values at
%                               frequency bin cutoffs with average of
%                               time-adjacent points (def: false)
%       'num_harmonics',[#]   - change the number of harmonics removed in
%                               deseasonalization process for daily data
%       'freq_bands',[names],[lims]  - set new endpoints for freq bands
%                                      default is [names] = {'< 5
%                                      days','5 - 30 days','30 - 365
%                                      days','> 1 year','all bands'} and
%                                      [lims] = [5;30;365;36500;0]
%       'multi_processor',[log],[#], - set whether to use multiple
%                                      processors and set number of [#]
%                                      processors to be used (def: true,
%                                      8)
%       'batch_processing',[log]     - sets job storage location
%                                      explicitly at
%                                      /tmp/kschwarzwald/[SLURMJOBID]
%                                      (def: false; only useful in batch
%                                      processing)
%       'save_log'            - exports log of attempted saves with
%                               success/failure status under
%                               ~/CalcLogs/Variability_CI_log_[starttime].txt
%       'subset_data',[type],[max1],[max2] - allows subsetting of data for
%                                            debugging/sampling process.
%                                            [type] = 'from1' will subset
%                                            lon from 1 to max1 and
%                                            lat from 1 to max2; [type] =
%                                            'rand' will subset [max1] lon
%                                            idxs and [max2] lat idxs
%                                            chosen at random based on the
%                                            [rand_seed] from above
%                                            (def:false)
%       'testing'             - adds '_TEST' to filename
%
%   Saving convention:
%   [filevar]_[freq]_[model]_[exp1]_[exp2]_[run]_StdDevsCI.mat
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   NOTE: THE CONFIDENCE INTERVALS ARE TAKEN FROM THE PURE RATIO
%
%   See also VARIABILITY, VARIABILITY_SEASONS
%   
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 12/12/2016 PARTIALLY IN BETA

%Load defaults
various_defaults = matfile('various_defaults.mat');

%% 1 Misc Setup
% Get timestamp (to identify logs)
format shortG
startTimestamp = clock;
start_msg = ['\nMeaningfulness and variability calculations run, \nstarted at ',...
                num2str(startTimestamp(4)),':',num2str(startTimestamp(5)),':',num2str(fix(startTimestamp(6))),...
                ' on ',num2str(startTimestamp(2)),'/',num2str(startTimestamp(3)),'/',num2str(startTimestamp(1)),...
                ',\nfor variables [',num2str(VarIndices),'] \nand models ',sprintf('%s ',modelArray{:}),'\nGood Luck! :) \n \n \n'];
fprintf(start_msg)

%Save command window output
diary_filename = [various_defaults.log_dir,'Variability_CI_out',num2str(startTimestamp(1)),num2str(startTimestamp(2)),num2str(startTimestamp(3)),num2str(startTimestamp(4))];
diary(diary_filename)

fullTime = tic;

%% 2 Set defaults
%Computing settings
batch_processing = false; 
multi_proc = true; 
num_proc = 8; %Apparently 8 is an optimized number

%Bootstrap settings
block_size = 365; 
seed = 65987436; 
nruns = 1200; %(block size = 400? boostrap num/nruns = 800?)

%Deseasonalizing and FFT settings
num_harmonics = 12; 
remove_locavg = false;

%Other defaults
expArray = {'rcp85','piControl'};
save_log = false; filename_add = [];
%pc_jobstorloc = strcat('/tmp/kschwarzwald/', getenv('SLURM_JOB_ID'));
pc_jobstorloc = strcat(various_defaults.tmp_dir,randi(10000));
subset_data = false;
replace_files = true;

%Define Frequency Bins - Period (days)
freq_band_setup = various_defaults.freq_band_setup;
freq_char_id = find(~cellfun(@isempty,strfind(freq_band_setup(:,1),'day')));
freq_bands = freq_band_setup{freq_char_id,3};
        freq_band_names = freq_band_setup{freq_char_id,2};
        freq_band_ID = freq_band_setup{freq_char_id,4};
num_bands = length(freq_bands);

%% 3 Set behavior of optional flags
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'experiments'}
                expArray = varargin{in_idx+1}; varargin{in_idx+1} = 0;
                if length(expArray)~=2
                    error('VariabilityCI:MissingExp',['Two experiments are needed for confidence intervals on ratios calculations, but ',num2str(length(expArray)),' were inputted.'])
                end
            case {'block_size'}
                block_size = varargin{in_idx+1};
            case {'nruns'}
                nruns = varargin{in_idx+1};
            case {'rand_seed'}
                seed = varargin{in_idx+1};
            case {'freq_bands'}
                freq_band_names = varargin{in_idx+1}; varargin{in_idx+1} = 0;
                freq_bands = varargin{in_idx+2}; varargin{in_idx+2} = 0;
                freq_band_ID = varargin{in_idx+3}; varargin{in_idx+3} = 0;
                num_bands = length(freq_bands);
            case {'multi_processor'}
                multi_proc = varargin{in_idx+1};
                if multi_proc == true
                    num_proc = varargin{in_idx+2};
                else
                	num_proc = 0;
                end
            case {'remove_locavg'}
                remove_locavg = varargin{in_idx+1};
            case {'batch_processing'}
                batch_processing = true;
            case {'num_harmonics'}
                num_harmonics = varargin{in_idx+1};
                if isa(varargin{in_idx}+1,'integer') == 0 || varargin{in_idx+1}<0
                    error('num_harmonics must be a positive integer')
                end
            case {'save_log'}
                save_log = true;
                complete_log = cell(length(modelArray)*length(VarIndices),1);
            case {'replace'}
                replace_files = varargin{in_idx+1};
            case {'subset_data'}
                subset_data = true;
                subset_type = varargin{in_idx+1};
                lon_lim = varargin{in_idx+2};
                lat_lim = varargin{in_idx+3};
            case {'testing'}
                filename_add = '_TEST_new';
        end
    end
end

                
%Set random number generator so that results are consistent between runs
rng(seed,'twister');

%% 4 Begin Parallel Computing
if multi_proc
    if isempty(gcp('nocreate'))
        % create a local cluster object
        pc = parcluster('local');
        
        %experimentlicitly set the JobStorageLocation to the temp directory that was
        %created in sbatch script (for batch processing only)
        if batch_processing
            pc.JobStorageLocation = pc_jobstorloc;
        end
        
        % start the matlabpool with num_workers workers
        parpool(pc, num_proc);
    end
end

%% 5 Compute CI
for var_idx = 1:length(VarIndices)
    %Get variable characteristics
    [~,filevar,freq,~,~,~,~,~] = var_chars(VarIndices(var_idx));
    for model_idx=1:length(modelArray)
        try
            %% 5.1 Model-Specific Graph Setup
            modelTime = tic;
            tic
            %Save variables needed across models and variables (to save
            %from clearing at end)
            initial_vars = who;
            
            %Isolate model
            model = modelArray{model_idx};
            fprintf(['\nCalculations for ',model,' beginning\n']);
            
            %Determine final filename
            [~,run,strtyr,endyr] = name_chars(model,expArray,filevar,freq,'use_convention',1,1,1);
            filenameS = [various_defaults.proc_data_dir,model,'/',filevar,freq,model,'_',expArray{1},'_',expArray{2},'_',run{1},'StdDevCI',filename_add,'.mat'];
            
            %% 5.2 Decide whether calculations or replacements are needed
            if replace_files %If replace_files, calculate and fully replace file
                calc = true; replace = true;
            else %If ~replace_files, find out whether the file exists, and if so, whether the wanted bands exist
                if exist(filenameS,'file')
                    load(filenameS);
                    bands_exist = cell(num_bands,1);
                    for band_test = 1:num_bands
                        bands_exist{band_test} = structfind(StdDevs,'Size',freq_bands(band_test,:));
                    end
                    if any(cellfun('isempty',bands_exist)) %If not all wanted bands exist, calculated, but add to file instead of replacing it
                        calc = true; replace = false;
                        if ~all(cellfun('isempty',bands_exist)) %Subset to only needed bands
                            [~,idx_tmp] = ismember(daily_freq_bands,cell2mat({StdDevs.Size}'),'rows'); freq_idx = find(~idx_tmp); clear idx_tmp
                            %Subset freq bands to calculate to just those
                            %not yet represented
                            freq_bands = freq_bands(freq_idx,:);
                            freq_band_names = freq_band_names(freq_idx);
                            freq_band_ID = freq_band_ID(freq_idx);
                            num_bands = size(freq_bands,1);
                        end
                    else %If wanted bands all exist, no need to do the thing
                        calc = false;
                    end
                else %If file doesn't exist, calculate and create file
                    calc = true; replace = true;
                end
            end
            
            %% 5.2 Execute Run if Desired
            if calc
                %% 5.2.1 Load Data, Initialize, and Get Data Characteristics
                %Get filenames
                filename1 = strcat(various_defaults.raw_data_dir,model,'/',filevar,freq,model,'_',expArray{2},'_',run{2},num2str(strtyr{2}),'-',num2str(endyr{2}),'_DATA.mat');
                filename2 = strcat(various_defaults.raw_data_dir,model,'/',filevar,freq,model,'_',expArray{1},'_',run{1},num2str(strtyr{1}),'-',num2str(endyr{1}),'_DATA.mat');
                Var1File = matfile(filename1); Var2File = matfile(filename2);
                
                %Load data, subsetting if desired
                if subset_data && strcmp(subset_type,'from1') %If subsetting in equally-spaced intervals (matfile doesn't support loading otherwise) (saves memory)
                    lon_idxs = (1:lon_lim);
                    lat_idxs = (1:lat_lim);
                    Raw_tmp{1} = Var1File.Raw(lon_idxs,lat_idxs,:);
                    Raw_tmp{2} = Var2File.Raw(lon_idxs,lat_idxs,:); clear Var2File
                    lat = Var1File.lat(lat_idxs,1);
                    lon = Var1File.lon(lon_idxs,1); clear Var1File
                else %If can't load from matfile
                    clear Var1File Var2File
                    load(filename1);
                    if subset_data && strcmp(subset_type,'rand') %if random subsetting THIS DOESN'T YET WORK - TO MAKE THIS WORK, HAVE TO USE LINEAR INDEXING
                        lon_idxs = randi(size(Raw,1),lon_lim,1);
                        lat_idxs = randi(size(Raw,2),lat_lim,1);
                    elseif subset_data && strcmp(subset_type,'region') %if regional subsetting
                        %Seems like a good thing to be able to quickly
                        %subset
                    elseif ~subset_data %If ~subset_data, just take entire arrays
                        lon_idxs = (1:size(Raw,1));
                        lat_idxs = (1:size(Raw,2));
                    end
                    Raw_tmp{1} = Raw(lon_idxs,lat_idxs,:); clear Raw
                    lat = lat(lat_idxs);
                    lon = lon(lon_idxs);
                    load(filename2);
                    Raw_tmp{2} = Raw(lon_idxs,lat_idxs,:); clear Raw
                    
                end
                %Get array sizes
                nlon=length(lon); nlat=length(lat);
                rge=size(Raw_tmp{1},3);

                Raw = Raw_tmp; clear Raw_tmp
                
                %Store start and end year, experiment, etc. characteristics
                %(since the new convention removes them from the filename)
                for experiment = 1:2
                    file_params(experiment).exp = expArray{experiment}; %#ok<AGROW>
                    file_params(experiment).strtyr = strtyr{experiment}; %#ok<AGROW>
                    file_params(experiment).endyr = endyr{experiment}; %#ok<AGROW>
                    file_params(experiment).run = run{experiment}; %#ok<AGROW>
                end
                
                %Initialize structure for memory efficiency
                StdDevsRatio_tmp(num_bands).Name = 'Initialize'; %#ok<AGROW>
                
                %Pre-allocate arrays
                run_idxs = zeros(nruns,rge);
                %Get indices of blocks to bootstrap through rng process
                parfor (l=1:nruns,num_proc)
                    %Get starting indices of each block in run
                    block_start_idxs = ceil(rge*rand(ceil(rge/block_size),1));
                    %Add indices for total of each block (block start idx +
                    %0:block_size-1)
                    blocks_idxs = bsxfun(@plus, block_start_idxs, 0:block_size-1)';
                    %Add those indices to the full bootstrap idx array
                    run_idxs(l,:) = blocks_idxs(1:rge);
                end
                clear blocks_idxs block_start_idxs
                
                toc
                intMsg = ['Setup for ',filevar,freq,model,' complete'];
                disp(intMsg)
                tic
                
                %% 5.2.2 Detrend and Deseasonalize By Longitude Bands
                %Calculate matrix for deseasonalization
                X0 = ones(rge,1);
                Xc = cos((2*pi/365)*(0:rge-1)'*(1:num_harmonics));
                Xs = sin((2*pi/365)*(0:rge-1)'*(1:num_harmonics));
                X =  [X0 Xc Xs];
                
                for experiment=1:2
                    %For each experiment, get time series
                    Raw_tmp = Raw{experiment};
                    
                    parfor (l=1:nlon,num_proc)
                        %Load and Detrend
                        Y = detrend(squeeze(Raw_tmp(l,:,:))')'; %MUST ' IN DETREND TO MAKE SURE DETRENDING HAPPENS BY RIGHT COLUMN (ACROSS TIME DIMENSION, NOT LATS)
                        %Deseasonalize
                        beta_hat = (X'*X)\(X'*Y');
                        Y = Y - beta_hat'*X';
                        %Store detrended, deseasonalized time series in
                        %original array
                        Raw_tmp(l,:,:) = Y;
                    end
                    Raw{experiment} = Raw_tmp;
                    clear yy_circ_tmp Raw_tmp Y
                end
                clear X0 Xc Xs X
                toc
                disp('Detrend, deseasonalization complete')
                tic
                
                %% 5.2.3 Bootstrapping
                %Get ratio of variability of each constructed bootstrap run
                band_stddev = cell(2,1);
                for experiment = 1:2
                    %Double the time series such that blocks can 'wrap around'
                    Raw_circ = cat(3,Raw{experiment},Raw{experiment}(:,:,1:block_size));
                    band_stddev_tmp = zeros(nruns,nlon,nlat,num_bands);
                    Ks = (0:(rge/2))./rge;
                    SFs = (1:num_harmonics)./365;
                    
                    %Put frequencies in order
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
                    end
                    
                    parfor (i=1:nlon,num_proc)  
                        Raw_circ_tmp = squeeze(Raw_circ(i,:,:));
                        for l=1:nruns
                            Y_c2 = rge.*abs(fft(squeeze(Raw_circ_tmp(:,run_idxs(l,:))),rge,2)./rge).^2; %#ok<PFBNS>
                            %Remove estimates at seasonal harmonics, replaces with local average
                            if remove_locavg
                                for x=1:length(SFs)
                                    [~,idx] = min(abs(SFs(x)-Ks));
                                    Y_c2(:,idx) = ( Y_c2(:,idx-1)+ Y_c2(:,idx+1))/2;
                                end
                                Y_c2(:,1)=( Y_c2(:,2)+ Y_c2(:,3))/2;
                            end

                            %Calculate and store band-separated standard
                            %deviations
                            for band = 1:num_bands
                                band_stddev_tmp(l,i,:,band) = ((sum(Y_c2(:,freq_bands_idxs(band,2):freq_bands_idxs(band,1)),2))/(rge/2)).^0.5; %#ok<PFBNS>
                            end
                        end
                    end
                    band_stddev{experiment} = band_stddev_tmp; clear band_stddev_tmp Raw_circ_tmp
                end
                
                %Calculate ratio of stddev of bootstrapped bands from above
                for band = 1:num_bands
                    StdDevsRatio_tmp(band).std_runs = squeeze(std(squeeze(band_stddev{1}(:,:,:,band))./squeeze(band_stddev{2}(:,:,:,band)),0,1));
                end
                clear band_stddev run_idxs
                toc
                progMsg = 'Bootstrapping complete';
                disp(progMsg);
                
                %% 5.2.4 Variability calc of original data
                tic
                %Get ratios of variability of the data (not bootstrap)
                for experiment = 1:2
                    %If variability data already exists, load it instead of
                    %recalculating it
                    filename1SD = [various_defaults.proc_data_dir,model,'/',filevar,freq,model,'_',expArray{experiment},'_',run{1},num2str(strtyr{experiment}),'-',num2str(endyr{experiment}),'_sqrtPower.mat'];
                    process_var = false; %Default, if file already exists. If tree below decides if not
                    if exist(filename1SD,'file')
                        test_vars = whos('-file',filename1SD);
                        if structfind(test_vars,'name','StdDevs')>0
                            test_vars_tmp = load(filename1SD,'StdDevs');
                            idx_size = zeros(num_bands,1);
                            for band = 1:num_bands
                                if ~isempty(structfind(test_vars_tmp.StdDevs,'Size',freq_bands(band,:)))
                                    idx_size(band) = structfind(test_vars_tmp.StdDevs,'Size',freq_bands(band,:));
                                else
                                    idx_size(band) = 0;
                                end
                            end
                            if any(idx_size==0)
                                process_var = true;
                            end
                        else %If in old file convention
                            process_var = true;
                        end
                    else
                        process_var = true;
                    end
                    
                    %If variability data doesn't yet exist, calculate it
                    if process_var 
                        disp(['Running Variability for ',expArray{experiment},' bands ',strjoin(freq_band_ID(~idx_size))])
                        Variability(VarIndices(var_idx),{model},'experiments',expArray(experiment),...
                            'daily_bands',freq_band_names(~idx_size),freq_bands(~idx_size,:),freq_band_ID(~idx_size),'replace',false);
                        test_vars_tmp = load(filename1SD,'StdDevs');
                    end
                    
                    for band = 1:num_bands
                        StdDevs_tmp(experiment,band).Data = test_vars_tmp.StdDevs(structfind(test_vars_tmp.StdDevs,'ID',freq_band_ID{band})).Data(lon_idxs,lat_idxs,:); %#ok<AGROW>
                    end
                    clear test_vars_tmp
                end
                
                clear period nfreqs Ks Raw_tmp Y_c2 SFs adj_freq_bands
                
                for band = 1:num_bands
                    StdDevsRatio_tmp(band).Name = freq_band_names{band};
                    StdDevsRatio_tmp(band).ID = freq_band_ID{band};
                    StdDevsRatio_tmp(band).Size = freq_bands(band);
                    StdDevsRatio_tmp(band).Data = StdDevs_tmp(1,band).Data./StdDevs_tmp(2,band).Data;
                end
                clear StdDevs_tmp
                toc
                disp(strcat('variability calcs complete'))
                
                %% 5.2.5 Determine unmeaningful pixels
                tic
                %Get indices of ratios not meaningfully different from 1
                for band = 1:num_bands
                    [dev_lon_tmp, dev_lat_tmp] = find(StdDevsRatio_tmp(band).Data - 2*StdDevsRatio_tmp(band).std_runs < 1 & StdDevsRatio_tmp(band).Data + 2*StdDevsRatio_tmp(band).std_runs > 1);
                    StdDevsRatio_tmp(band).lon_ciidx = dev_lon_tmp; clear dev_lon_tmp
                    StdDevsRatio_tmp(band).lat_ciidx = dev_lat_tmp; clear dev_lat_tmp
                end
                disp('meaningful pixels determined')
                
                %% 5.2.6 Save output
                %Add calculated bands to existing ones, if desired
                if ~replace
                    StdDevsRatio = [StdDevsRatio';StdDevsRatio_tmp']'; clear StdDevsRatio_tmp;
                else
                    StdDevsRatio = StdDevsRatio_tmp; clear StdDevsRatio_tmp
                end
                
                save(filenameS,'-v7.3','StdDevsRatio','file_params','lat','lon');
                
                toc
                disp(['*****CI for ',model,' ',filevar,freq,' Complete*****']);
                clear run_idxs lat lon
                
                if save_log
                    save_log_msg = [filevar,freq,model,' Complete'];
                    complete_log{(var_idx-1)*length(modelArray)+model_idx} = save_log_msg;
                end
                
            else %If file already exists
                disp(['Files for ',model,' ',expArray{1},' / ',expArray{2},' ',filevarFN,freq,' already exist, were not replaced']);
                if save_log
                    save_log_msg = [filevar,freq,model,' already exists'];
                    complete_log{(var_idx-1)*length(modelArray)+model_idx} = save_log_msg;
                end
            end
            
        catch ME %If errors
            disp(ME);
            if save_log
                save_log_wrn = ['****',filevar,freq,model,' Not Complete****'];
                complete_log{(var_idx-1)*length(modelArray)+model_idx} = save_log_wrn;
            end
            try  %#ok<TRYNC>
                for a = 1:length(ME.stack)
                    disp(ME.stack(a))
                end
            end
            try disp(ME.stack.line); end %#ok<TRYNC>
            try disp(ME.cause); end %#ok<TRYNC>
        end
        modelEnd = toc(modelTime);
        disp(['total time for ',model,' ',filevar,freq,': ',num2str(modelEnd)])
        clearvars('-except',initial_vars{:});
    end
    
    %Export log as table
    if save_log
        export_log = cell2table(complete_log,'VariableNames',{'CI_attempted_executions'});
        log_filename = [various_defaults.log_dir,'Variability_CI_log',num2str(startTimestamp(1)),'-',num2str(startTimestamp(2)),'-',num2str(startTimestamp(3)),'-',num2str(startTimestamp(4)),'-',num2str(startTimestamp(5)),'-',num2str(startTimestamp(6)),'.txt'];
        writetable(export_log,log_filename)
    end
    disp(['********CI  for ',filevar,freq,' for All Models Complete********']);
end

%% 6 Finish diary, clock total project
fullEnd = toc(fullTime);
disp(['total time for all models, variables: ',num2str(fullEnd)])
endTimestamp = clock;
disp(endTimestamp)
diary off

end
