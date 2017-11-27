function SpectralRatios(VarIndices,modelArray,varargin)
% SPECTRALRATIOS   calculate spectral ratios.
%   SPECTRALRATIOS(VarIndices,modelArray) calculates the simple ratio of
%   variability split up by frequency band between two experiments
%   (default: RCP8.5 / piControl) for variables given by [VarIndices] and
%   models given by [modelArray]. The frequency bands to be processed are
%   all those common between the loaded files from both experiments and are
%   taken from the struct StdDev from '_sqrtPower.mat' files. The resultant
%   ratios are saved in a struct [StdDevRatios] with a length equal to the
%   number of frequency band sizes common to both StdDev structs, and the
%   corresponding .Name, .ID, and .Size fields as the StdDev struct from
%   the first experiment. Files are saved in [proc_data_dir]/[model]/
%   for full-year data and [proc_data_dir]/[model]/[season_dir]/
%   for seasonal data (see below).
%
%   Function will attempt to use variability calculated for two runs of the
%   same length, will throw warning if run lengths are different.
%
%   SPECTRALRATIOS(...,'[flag]',[params],...) modify program run as below:
%       'experiments',[cell]- manually set experiments to process, first
%                             will be divided by second (default
%                             {'rcp85','piControl'})
%       'replace',[log]     - set whether to replace existing files (def:
%                             true)
%       'seasons',[seasons],([various_defaults.season_dir])
%                           - process seasonal data for the seasons in the
%                             array [season] (i.e. {'JJA','DJF'}).
%                             Optionally change the season folder (in
%                             [proc_data_dir]/[model]) from which to
%                             get and save files by also setting
%                             [various_defaults.season_dir] = '/[folder name]/'
%       'break_on_length_mismatch',[log] - set whether to throw an error if
%                                          length of time between the two
%                                          experiments is not equal (if
%                                          false, then a warning is thrown
%                                          but the spectral ratio file is
%                                          still created), def: true. 
%
%
%   Saving convention:
%   [filevar]_[freq]_[model]_[exp1]_[exp2]_[run]_SpectralRatios(_[seas]).mat
%
%   Sample output (for each frequency band [band]):
%   StdDevRatios(band).Name = [char] (freq band description, i.e. '< 5
%                                     days')
%   StdDevRatios(band).ID = [char] (Identifier, i.e. 'HF','LF',...)
%   StdDevRatios(band).Size = [2 x 1 array] (freq band size, i.e. [0,5],...)
%   StdDevRatios(band).Data = [nlon x nlat array] (Variability Ratios)
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 10/14/2016

replace_files = true;
expArray = {'rcp85','piControl'};
seasArray = {'all'};
various_defaults = matfile('various_defaults.mat');
length_break = true;
cust_timeframe = [];
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'replace'}
                replace_files = varargin{in_idx+1};
            case {'experiments'}
                expArray = varargin{in_idx+1};
                varargin{in_idx+1} = 0; %(To get around switch tripping over the cell)
            case {'seasons'}
                seasArray = varargin{in_idx+1}; varargin{in_idx+1}=0;
                if length(varargin)>in_idx+1
                    if strcmp(varargin{in_idx+2}(1),'/')
                        various_defaults.season_dir = varargin{in_idx+2};
                    end
                end
            case {'timeframe'}
                cust_timeframe = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'break_on_length_mismatch'}
                length_break = varargin{in_idx+1};
        end
    end
end


for i = 1:length(VarIndices)
    %Get variable characteristics
    [~,filevarFN,freq,~,~,~,~,~,~] = var_chars(VarIndices(i));
    for j=1:length(modelArray)
        for seas = 1:length(seasArray)
            intermod_vars = who;
            try
                %Set filename to include season, if subset by season
                if strcmp(seasArray{seas},'all')
                    seas_filen = []; seas_folder = '/';
                    [~,run,strtyr,endyr] = name_chars(modelArray{j},expArray,filevarFN,freq);
                else seas_filen = ['_',seasArray{seas}]; seas_folder = various_defaults.season_dir;
                    [~,run,strtyr,endyr] = name_chars(modelArray{j},expArray,filevarFN,freq,'season',seasArray{seas},various_defaults.season_dir);
                end
                
                if isempty(endyr{1}) || isempty(endyr{2})
                    error('SpectralRatios:MissingFile','Missing file (see name_chars warning above), no file created.')
                end
                
                if length_break && (endyr{2}-strtyr{2} ~= endyr{1}-strtyr{1})
                    error('SpectralRatios:LengthMismatch','No set of two runs from the inputted experiments were found that were the same length, no file created.') 
                end
                
                %Get file names/characteristics of files of inputted var/exp/model
                warning('off','name_chars:MultConvs'); %Turn off multiple convention warning (it's annoying)
                
                %Store start and end year, experiment, etc. characteristics
                %(since the new convention removes them from the filename)
                for experiment = 1:2
                    file_params(experiment).exp = expArray{experiment}; %#ok<AGROW>
                    file_params(experiment).strtyr = strtyr{experiment}; %#ok<AGROW>
                    file_params(experiment).endyr = endyr{experiment}; %#ok<AGROW>
                    file_params(experiment).run = run{experiment}; %#ok<AGROW>
                end
                
                %% Process Data
                %Set filename to save ratios in
                filenameS = [various_defaults.proc_data_dir,modelArray{j},seas_folder,filevarFN,freq,modelArray{j},'_',expArray{1},'_',expArray{2},'_',run{1},'SpectralRatios',seas_filen,'.mat'];
                %Calculate ratios, if necessary
                if replace_files || (~replace_files && ~exist(filenameS,'file'))
                    
                    %Load Standard Deviations for first experiment
                    StdDevs1 = load_stddevs(VarIndices(i),modelArray{j},expArray{1},seasArray{seas},'start_year',strtyr{1},'various_defaults.season_dir',various_defaults.season_dir);
                    
                    %Load Standard Deviations for second experiment
                    StdDevs2 = load_stddevs(VarIndices(i),modelArray{j},expArray{2},seasArray{seas},'start_year',strtyr{2},'various_defaults.season_dir',various_defaults.season_dir);
                    
                    %Load lat / lon values
                    [~,lat,lon] = load_stddevs(VarIndices(i),modelArray{j},expArray{1},seasArray{seas},'start_year',strtyr{1},'various_defaults.season_dir',various_defaults.season_dir); %#ok<ASGLU>
                    
                    if ~exist('StdDevs1','var') || ~exist('StdDevs2','var')
                        error('SpecRat:NoData','Required data not found, exiting')
                    end
                    
                    %Get size contents of each StdDevs struct
                    sizes1 = cell2mat({StdDevs1(1:length(StdDevs1)).Size}');
                    sizes2 = cell2mat({StdDevs2(1:length(StdDevs2)).Size}');
                    %Get ID contents of each StdDevs struct
                    ids1 = {StdDevs1(1:length(StdDevs1)).ID}';
                    ids2 = {StdDevs2(1:length(StdDevs2)).ID}';
                    %Get band sizes common to both
                    [size_intersect] = intersect(sizes1,sizes2,'rows');
                    [ids_intersect] = intersect(ids1,ids2);
                    if isempty(size_intersect)
                        error('SpecRat:NoMatchingBands','The two StdDevs structs do not have any matching length frequency bands, no ratios will be calculated')
                    end
                    %As default, check against periods per band
                    intersects = mat2cell(size_intersect,ones(size(size_intersect,1),1),size(size_intersect,2));
                    field_match = 'Size';
                    %If periods not specific enough (if, say, several
                    %calculations contained in struct over same frequency
                    %bands), use IDs if possible or throw warning
                    if length(size_intersect) ~= length(StdDevs1) || length(size_intersect) ~= length(StdDevs2)
                        if length(ids_intersect) == length(StdDevs1) && length(ids_intersect) == length(StdDevs2)
                            intersects = ids_intersect; field_match = 'ID';
                        else
                            warning('SpecRat:UnequalBands','The two StdDevs structs do not have data for all the same length frequency bands, ratios will only be calculated for those in common')
                        end    
                    end
                    
                    %Calculate ratios for band sizes common to both
                    for band = 1:length(intersects)
                        %find indices of each frequency band size to process
                        idx1 = structfind(StdDevs1,field_match,intersects{band});
                        idx2 = structfind(StdDevs2,field_match,intersects{band});
                        %Calculate spectral ratios
                        StdDevsRatio(band).Data = StdDevs1(idx1).Data./StdDevs2(idx2).Data; %#ok<AGROW>
                        %Calculate confidence intervals
                        if isfield(StdDevs1,'std_boot') && isfield(StdDevs2,'std_boot')
                            [dev_lon_tmp, dev_lat_tmp] = find(StdDevsRatio(band).Data - 2*StdDevsRatio(band).std_runs < 1 & StdDevsRatio(band).Data + 2*StdDevsRatio(band).std_runs > 1);
                            StdDevsRatio(band).lon_ciidx = dev_lon_tmp; clear dev_lon_tmp
                            StdDevsRatio(band).lat_ciidx = dev_lat_tmp; clear dev_lat_tmp
                        end
                        %Save name, id, and size
                        StdDevsRatio(band).Name = StdDevs1(idx1).Name; %#ok<AGROW>
                        StdDevsRatio(band).ID = StdDevs1(idx1).ID; %#ok<AGROW>
                        StdDevsRatio(band).Size = StdDevs1(idx1).Size; %#ok<AGROW>
                    end
                    
                    %Save Spectral Ratios
                    save(filenameS,'-v7.3','lat','lon','file_params','StdDevsRatio');
                    clear StdDevsRatio lat lon
                    disp([filenameS,' saved.'])
                    endmsg = ['Spectral ratios for ',modelArray{j},' ',filevarFN,freq,' ',seas_filen,' saved.'];
                    disp(endmsg)
                else
                    disp([filenameS,' already exists, was not replaced']);
                end
                clearvars('-except',intermod_vars{:});
            catch ME
                disp(ME.message); %disp(ME.stack);
                warningMsg = ['Spectral ratios not calculated for ',modelArray{j},' ',filevarFN,freq];
                warning(warningMsg)
                clearvars('-except',intermod_vars{:});
            end
        end
    end
end

end