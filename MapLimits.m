function MapLimits(VarIndices,modelArray,varargin)
% MAPLIMITS     Get intermodel unweighted means and standard deviations for
%               colorbar limiting purposes
%
%   MAPLIMITS(VarIndices,modelArray) calculates the global mean, standard
%   deviation, median, and median absolute deviation of each calculation
%   contained in the files _sqrtPower and _LocalMeans for each variable in
%   [VarIndices] and model in [modelArray]. Then, the intermodel mean of
%   each of those global mean values is calculated, outputted into a struct
%   [IntermodelStats] or [Intermodel_Ratio_Stats], and saved as detailed
%   below. This process is repeated for (by default) RCP8.5, and piControl
%   data and for the ratio RCP8.5/piControl of each dataset, intermodel
%   statistics of which are all saved in separate files [_IntermodelStats]
%   for the former and [_IntermodelRatioStats] for the latter. The structs
%   contain the fields [ID] (variable identifier, i.e. 'HF'), [mean] (means
%   of all global means), [std] (means of all global standard deviations),
%   [med] (means of all global medians), [mad] (means of all global median
%   absolute deviations), and [models] (the models from which the
%   intermodel statistic is taken). If a model doesn't have a certain
%   variable, a warning is thrown and that model is ignored for the
%   purposes of that variable's calculations, but program run continues.
%
%   NOTE: THESE INTERMODEL STATS ARE *NOT* SUMMARY STATISTICS - THEY ARE
%   SOLELY DESIGNED TO ALLOW FOR AN EASY WAY TO CREATE COLORMAPS THAT ARE
%   CONSISTENT ACROSS MODELS. As such, these statistics are not weighted by
%   area, wholly ignore NaNs, and make no complex assumptions about
%   underlying distributions. They are simply the mean of the means of the
%   dataset values across models and the mean of the standard deviations of
%   these dataset values across models.
%
%   Function requirements (on path): var_chars, name_chars, load_stddevs,
%   load_means, wraptext (from MathWorks community)
%
%   Sample output of MAPLIMITS(20,{'ACCESS1-3','BCC-CSM1-M','CanESM2',...})
%       1) One file [_IntermodelStats] each for RCP8.5 and piControl
%       containing: 
%         IntermodelStats.ID        = {'HF','MF','LF','XF','mean',...}
%         IntermodelStats.mean      = {[#] ,[#] ,[#] ,[#] ,[#],...}
%         IntermodelStats.std       = {[#] ,[#] ,[#] ,[#] ,[#],...}
%         IntermodelStats.med       = {[#] ,[#] ,[#] ,[#] ,[#],...}
%         IntermodelStats.mad       = {[#] ,[#] ,[#] ,[#] ,[#],...}
%         IntermodelStats.models    = {[ar],[ar],[ar],[ar],[ar],...}
%       2) One file [_IntermodelRatioStats], using mean, std from each
%       dataset in RCP8.5 divided by the corresponding in piControl,
%       containing:
%         IntermodelRatioStats.ID   = {'HF','MF','LF','XF','mean',...}
%         IntermodelRatioStats.mean = {[#] ,[#] ,[#] ,[#] ,[#],...}
%         IntermodelRatioStats.std  = {[#] ,[#] ,[#] ,[#] ,[#],...}
%         IntermodelRatioStats.med  = {[#] ,[#] ,[#] ,[#] ,[#],...}
%         IntermodelRatioStats.mad  = {[#] ,[#] ,[#] ,[#] ,[#],...}
%         IntermodelRatioStats.models={[ar],[ar],[ar],[ar],[ar],...}
%   for [ar] = {'ACCESS1-3','CanESM2',...}
%
%   MAPLIMITS(...,'[flag]',[params],...) modify program run as below:
%       'experiments',[array]       - manually set experiment(s) to
%                                     process. If two experiments are
%                                     chosen, then intermodel stats are 
%                                     calculated for the ratio of the first
%                                     to the second experiment as well
%                                     (def. {'rcp85','piControl'})
%       'season',[char]             - calculates instead for seasonal data
%                                     (i.e. 'JJA')
%       'season_folder',[char]      - changes folder in
%                                     [proc_data_dir]/[model] in which
%                                     seasonal data is looked for (def:
%                                     [season_dir])
%
%   Saving convention (in directory [proc_data_dir]/IntermodelStats/):
%   [filevar]_[freq]_[exp]_IntermodelStats(_[season])
%   [filevar]_[freq]_[exp1]_[exp2]_IntermodelRatioStats(_[season])
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also: MAPS_DIAGNOSTICS (for applications)
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 11/26/2017
%
%   SUGGESTION - SIMPLIFY CODE BY CHANGING ALL OF THIS TO A LOOP OVER
%   OPERATIONS (i.e. {'nanmean','nanstd','nanmedian','mad'}, for evals, but
%   then using dynamic struct names, IntermodelStats.([str]) to save). THIS
%   WOULD ALLOW FOR MUCH EASIER ADDING OF FUNCTIONALITIES IN THE FUTURE

%% 1 Intro
disp(wraptext(['Intermodel means and stds calculations beginning for models ',strjoin(modelArray)]))

%% 2 Set Defaults / Deal with Function Flags
various_defaults = matfile('various_defaults.mat');
expArray = {'rcp85','piControl'};
process_byseason = false;
folder_season = various_defaults.season_dir;

%Set optional flag behavior
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'experiments'}
                expArray = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'season'}
                process_byseason = true;
                season = varargin{in_idx+1};
            case {'season_folder'}
                folder_season = varargin{in_idx+1};
        end
    end
end

num_exps = length(expArray);
num_models = length(modelArray);

%% 3 Calculate
for i = 1:length(VarIndices)
    [~,filevarFN,freq] = var_chars(VarIndices(i));
    
    %% 3.1 Determine total number of variables to calculate stats of
    IDs = cell(length(modelArray),1);
    %Get frequencies present in the variables, models looked at
    for j = 1:length(modelArray)
        if ~process_byseason
            season_load = 'all';
            [~,~,strtyr,~,~] = name_chars(modelArray{j},expArray,filevarFN,freq);
        else
            season_load = season;
            folder = folder_season; dir_seas = [various_defaults.proc_data_dir,modelArray{j},folder];
            [~,~,strtyr,~,~] = name_chars(modelArray{j},expArray,filevarFN,freq,'directory',dir_seas,'file_end',['*Power_',season,'.mat']);
        end
        
        IDs{j} = load_stddevs(VarIndices(i),modelArray{j},expArray{1},season_load,'start_year',strtyr{1},'season_folder',folder_season,'property','ID');
    end
    
    %Get full list of stored frequency bands/calculations (without repeats)
    IDs = unique([IDs{:}]);
    num_freqs = length(IDs);
    
    %Pre-allocate final output structs, based on those IDs, adding the mean
    %values as well.
    IntermodelStats_tmp = struct();
    Intermodel_Ratio_Stats = struct('ID',{IDs,'mean'});
    
    %Pre-allocate for collecting data from different experiments
    StdDevs = cell(num_exps,1);
    Means = cell(num_exps,1);
    
    %Pre-allocate intermediate arrays to calculate
    mean_means = zeros(num_models,num_exps); mean_med = zeros(num_models,num_exps);
    mean_stds = zeros(num_models,num_exps); mean_mad = zeros(num_models,num_exps);
    std_means = zeros(num_freqs,num_models,num_exps); std_med = zeros(num_freqs,num_models,num_exps);
    std_stds = zeros(num_freqs,num_models,num_exps); std_mad = zeros(num_freqs,num_models,num_exps);
    
    mean_ratio_means = zeros(num_models,1); mean_ratio_med = zeros(num_models,1);
    mean_ratio_stds = zeros(num_models,1); mean_ratio_mad = zeros(num_models,1);
    std_ratio_means = zeros(num_freqs,num_models); std_ratio_med = zeros(num_freqs,num_models);
    std_ratio_stds = zeros(num_freqs,num_models); std_ratio_mad = zeros(num_freqs,num_models);

    %% 3.2 Calculate global mean, standard deviation of each model's datasets
    for j = 1:length(modelArray)
        if ~process_byseason
            season_load = 'all';
            [~,~,strtyr,~,~] = name_chars(modelArray{j},expArray,filevarFN,freq);
        else
            season_load = season;
            folder = folder_season; dir_seas = [various_defaults.proc_data_dir,modelArray{j},folder];
            [~,~,strtyr,~,~] = name_chars(modelArray{j},expArray,filevarFN,freq,'directory',dir_seas,'file_end',['*Power_',season,'.mat']);
        end
        
        for experiment = 1:num_exps
            %Load required data for first variable
            StdDevs{experiment} = load_stddevs(VarIndices(i),modelArray{j},expArray{experiment},season_load,'start_year',strtyr{experiment},'season_folder',folder_season);
            Means{experiment} = load_means(VarIndices(i),modelArray{j},expArray{experiment},season_load,'start_year',strtyr{experiment},'season_folder',folder_season);
            
            %Get means, std devs of the data contained in Means
            mean_means(j,experiment) = nanmean(Means{experiment}(:));
            mean_stds(j,experiment) = nanstd(Means{experiment}(:));
            
            %Get median, median absolute deviation of the data contained in
            %Means
            mean_med(j,experiment) = nanmedian(Means{experiment}(:));
            mean_mad(j,experiment) = mad(Means{experiment}(:));
            
            for band = 1:num_freqs
                %Get relevant index in the IDs cell
                band_idx = find(cellfun(@(x) strcmp(IDs{band},x),{StdDevs{experiment}(:).ID}));
                if ~isempty(band_idx)
                    %Get means and stds of the frequency band variability
                    %contribution values
                    std_means(band,j,experiment) = nanmean(StdDevs{experiment}(band_idx).Data(:));
                    std_stds(band,j,experiment) = nanstd(StdDevs{experiment}(band_idx).Data(:));
                    
                    %Get median, median absolute deviation of the data
                    %contained in StdDevs
                    std_med(band,j,experiment) = nanmedian(StdDevs{experiment}(band_idx).Data(:));
                    std_mad(band,j,experiment) = mad(StdDevs{experiment}(band_idx).Data(:));
                end
            end
        end
        
        %If two experiments chosen, get stats for ratios, too
        if num_exps==2
            %Get means, std devs of the data contained in the ratio of Means
            Means_Ratio = Means{1}./Means{2};
            mean_ratio_means(j) = nanmean(Means_Ratio(:));
            mean_ratio_stds(j) = nanstd(Means_Ratio(:));
            
            %Get medians, median absolute deviations of the data contained
            %in the ratio of Means
            mean_ratio_med(j) = nanmedian(Means_Ratio(:));
            mean_ratio_mad(j) = mad(Means_Ratio(:));
            
            for band = 1:num_freqs;
                band_idx1 = find(cellfun(@(x) strcmp(IDs{band},x),{StdDevs{1}(:).ID}));
                band_idx2 = find(cellfun(@(x) strcmp(IDs{band},x),{StdDevs{2}(:).ID}));
                if ~isempty(band_idx1) && ~isempty(band_idx2)
                    StdDevs_Ratio = StdDevs{1}(band_idx1).Data./StdDevs{2}(band_idx2).Data;
                    %Get means, std devs of the data contained in the ratio
                    %of StdDevs
                    std_ratio_means(band,j) = nanmean(StdDevs_Ratio(:));
                    std_ratio_stds(band,j) = nanstd(StdDevs_Ratio(:));
                    
                    %Get medians, median absolute deviations of the data
                    %contained in the ratio of StdDevs
                    std_ratio_med(band,j) = nanmedian(StdDevs_Ratio(:));
                    std_ratio_mad(band,j) = mad(StdDevs_Ratio(:));
                end
            end
        end
    end
    
    %% 3.3 Calculate intermodel mean and standard devation of each global mean/std 
    for experiment = 1:num_exps
        initial_vars = who; %Get stock, to delete all variables created in this loop at the end, just for safety
        models = cell(num_freqs+1,1);
        for band = 1:num_freqs
            std_means_tmp = squeeze(std_means(band,:,experiment));
            std_stds_tmp = squeeze(std_stds(band,:,experiment));
            std_meds_tmp = squeeze(std_med(band,:,experiment));
            std_mads_tmp = squeeze(std_mad(band,:,experiment));
            %Calculate intermodel stats for struct data
            IntermodelStats_tmp(experiment,band).mean = mean(std_means_tmp(std_means_tmp ~= 0));
            IntermodelStats_tmp(experiment,band).std = mean(std_stds_tmp(std_stds_tmp ~= 0));
            IntermodelStats_tmp(experiment,band).med = mean(std_meds_tmp(std_meds_tmp ~= 0));
            IntermodelStats_tmp(experiment,band).mad = mean(std_mads_tmp(std_mads_tmp ~= 0));
            if ~isempty(std_means_tmp(std_stds_tmp==0))
                warning('MapLimits:missing_models',['Intermodel calculations for ',expArray{experiment},' ',IDs{band},' do not take into account ',strjoin(modelArray(std_stds_tmp == 0)),'.'])
            end
            models{band} = modelArray(std_stds_tmp ~= 0);
        end
        %Calculate intermodel stats for mean data
        mean_means_tmp = squeeze(mean_means(:,experiment));
        mean_stds_tmp = squeeze(mean_stds(:,experiment));
        mean_meds_tmp = squeeze(mean_med(:,experiment));
        mean_mads_tmp = squeeze(mean_mad(:,experiment));
        IntermodelStats_tmp(experiment,num_freqs+1).mean = mean(mean_means_tmp(mean_means_tmp ~= 0));
        IntermodelStats_tmp(experiment,num_freqs+1).std = mean(mean_stds_tmp(mean_stds_tmp ~= 0));
        IntermodelStats_tmp(experiment,num_freqs+1).med = mean(mean_meds_tmp(mean_meds_tmp ~= 0));
        IntermodelStats_tmp(experiment,num_freqs+1).mad = mean(mean_mads_tmp(mean_mads_tmp ~= 0 ));
        if ~isempty(mean_stds_tmp(mean_stds_tmp==0))
            warning('MapLimits:missing_models',['Intermodel calculations for ',expArray{experiment},' means do not take into account ',strjoin(modelArray(mean_means_tmp == 0)),'.'])
        end
        models{num_freqs+1} = modelArray(mean_means_tmp ~= 0);
        
        %Get out intermodel stats for this experiment
        IntermodelStats = squeeze(IntermodelStats_tmp(experiment,:));
        %Add frequency/data type idenitifiers
        IDs_and_mean = [IDs,{'mean'}];
        [IntermodelStats(:).ID] = IDs_and_mean{:};
        [IntermodelStats(:).models] = models{:};
        
        if ~exist([various_defaults.proc_data_dir,'/IntermodelStats/'],'file')
            mkdir([various_defaults.proc_data_dir,'/IntermodelStats/'])
            warning([various_defaults.proc_data_dir,'/IntermodelStats/ created.'])
        end
        
        filenameS = [various_defaults.proc_data_dir,'/IntermodelStats/',filevarFN,freq,expArray{experiment},'_IntermodelStats'];
        if process_byseason
            filenameS = [filenameS,'_',season]; %#ok<AGROW>
        end
        save(filenameS,'-v7.3','IntermodelStats')
        disp([filenameS,' saved!'])
        clearvars('-except',initial_vars{:});
    end
    
    %If two experiments, do the same for the ratio
    if num_exps==2
        models = cell(num_freqs+1,1);
        for band = 1:num_freqs
            std_ratio_means_tmp = squeeze(std_ratio_means(band,:));
            std_ratio_stds_tmp = squeeze(std_ratio_stds(band,:));
            std_ratio_mads_tmp = squeeze(std_ratio_mad(band,:));
            std_ratio_meds_tmp = squeeze(std_ratio_med(band,:));
            %Calculate intermodel stats for struct data
            Intermodel_Ratio_Stats(band).mean = mean(std_ratio_means_tmp(std_ratio_means_tmp ~= 0));
            Intermodel_Ratio_Stats(band).std = mean(std_ratio_stds_tmp(std_ratio_stds_tmp ~= 0));
            Intermodel_Ratio_Stats(band).mad = mean(std_ratio_mads_tmp(std_ratio_mads_tmp ~= 0));
            Intermodel_Ratio_Stats(band).med = mean(std_ratio_meds_tmp(std_ratio_meds_tmp ~= 0));
            if ~isempty(std_ratio_means_tmp(std_ratio_stds_tmp==0))
                warning('MapLimits:missing_models',['Intermodel calculations for ',strjoin(expArray(:)),' ratio ',IDs{band},' do not take into account ',strjoin(modelArray(std_ratio_stds_tmp == 0)),'.'])
            end
            models{band} = modelArray(std_ratio_stds_tmp ~= 0);
        end
        %Calculate intermodel stats for mean data
        Intermodel_Ratio_Stats(num_freqs+1).mean = mean(mean_ratio_means(mean_ratio_means ~= 0));
        Intermodel_Ratio_Stats(num_freqs+1).std = mean(mean_ratio_stds(mean_ratio_stds ~= 0));
        Intermodel_Ratio_Stats(num_freqs+1).med = mean(mean_ratio_med(mean_ratio_med ~= 0));
        Intermodel_Ratio_Stats(num_freqs+1).mad = mean(mean_ratio_mad(mean_ratio_mad ~= 0));
        if ~isempty(mean_ratio_stds(mean_ratio_stds==0))
            warning('MapLimits:missing_models',['Intermodel calculations for ',strjoin(expArray(:)),' ratio means do not take into account ',strjoin(modelArray(mean_ratio_means == 0)),'.'])
        end
        models{num_freqs+1} = modelArray(mean_ratio_means ~= 0);
        
        %Add frequency/data type idenitifiers
        IDs_and_mean = [IDs,{'mean'}];
        [Intermodel_Ratio_Stats(:).ID] = IDs_and_mean{:};
        [Intermodel_Ratio_Stats(:).models] = models{:};
        
        filenameS = [various_defaults.proc_data_dir,'IntermodelStats/',filevarFN,freq,expArray{1},'_',expArray{2},'_IntermodelRatioStats'];
        if process_byseason
            filenameS = [filenameS,'_',season]; %#ok<AGROW>
        end
        save(filenameS,'-v7.3','Intermodel_Ratio_Stats')
        disp([filenameS,' saved!'])
    end
    
end

end