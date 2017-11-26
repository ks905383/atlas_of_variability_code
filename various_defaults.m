% VARIOUS_DEFAULTS  set conventions across the file ecosystem
%
%   Conventions/variables included in the file various_defaults.mat:
%       - code_dir          - directory in which code is stored.
%       - raw_data_dir      - parent directory in which raw (netCDF) data
%                             is stored (data is stored in subfolders
%                             marked by model name).
%       - proc_data_dir     - parent directory in which processed data is
%                             stored (data is stored in subfolders marked
%                             by model name). Can be the same as
%                             raw_data_dir above.
%       - season_dir        - subfolder of [proc_data_dir]/[model]/
%                             containing seasonal data. It is recommended
%                             that it is kept separate from the processed
%                             data directory (by having it not set to []).
%       - figure_dir        - directory under which output figures are
%                             saved, stored by 'domain' in subdirectories
%                             (i.e. 'Precipitation').
%       - gridweights_dir   - directory that contains the area weights for
%                             climate models. Can be the same as any of the
%                             above.
%                            
%
%       - expArray_disp     - (N x 2) cell array containing the 'raw'/save
%                             name of used experiments in the first column
%                             (i.e. 'i1400') and the display name of those
%                             experiments (i.e. '1400 ppm') in the second
%                             column. Sample code to extract the display
%                             name using an [expArray] cell array: exp_name = various_defaults.expArray_disp(find(cellfun(@(x) strcmp(expArray{1},x),various_defaults.expArray_disp(:,1))),2);      
%       - modelArray_disp   - (N x 2) cell array containing the 'raw'/save
%                             name of used experiments in the first column
%                             (i.e. 'CSIRO-Mk3-6-0') and the display/actual
%                             name of those models (i.e. CSIRO Mk3.6.0) in
%                             the second column. Only models that have a
%                             discrepancy in the two are used. Sample code
%                             to extract the display name using a
%                             [modelArray] cell array: model_name = various_defaults.modelArray_disp(find(cellfun(@(x) strcmp(modelArray{1},x),various_defaults.modelArray_disp(:,1))),2);
%       
%       - freq_band_setup   - cell array containing the default frequency
%                             bands processed in the program Variability.m.
%                             See the relevant code section below for
%                             format. Variability.m will by default process
%                             the frequency bands listed here.
%
%       - variability_ids   - the IDs of statistics of variability (versus
%                             the mean). This is used in Maps_general to 
%                             determine some graph titling conventions. 
%
%   Note: this program is part of the /project/moyer/ climate data file
%   ecosystem. 
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 11/15/2017 

clear
%% Default locations 
%Basic code directory
code_dir = '~/Documents/GitHub/atlas_of_variability_code/';

%Raw Data directory
raw_data_dir = '~/Documents/climate_data/';

%Processed Data directory
proc_data_dir = '~/Documents/climate_data/';

%Seasonal Data directory (in each [proc_data_folder]/[model]/)
season_dir = '/seasonal_data/';

%Figure save directory
figure_dir = '~/Documents/climate_data/';

%Gridweights directory
gridweights_dir = '~/Documents/climate_data/';

%Temporary storage directory (for multi-processor computing)
tmp_dir = '~/Documents/climate_data/';

%Log storage directory (for logs of code runs)
log_dir = '~/Documents/climate_data/';

%% Default expressions
%Display names of experiment runs (column 1 = 'raw'/save name, column 2 =
%display name)
expArray_disp = {'rcp85','RCP8.5';...
                 'piControl','PI';...
                 'historical','historical';...
                 'base','present';...
                 'i1400','1400 ppm';...
                 '1870control','PI';...
                 'rcp85REG','RCP8.5';...
                 'piControlREG','PI';...
                 'i1400REG','1400 ppm';...
                 '1870controlREG','PI';...
                 'historical1','historical, first half';...
                 'historical2','historical, second half';...
                 'abrupt4x','abrupt 4x CO2';...
                 'abrupt2x','abrupt 2x CO2';...
                 'abrupt8x','abrupt 8x CO2';...
                 'abrupt6x','abrupt 6x CO2';...
                 'abrupt16x','abrupt 16x CO2';...
                 'abrupt32x','abrupt 32x CO2';...
                 'control','PI'};

%Display names of model runs (column 1 = 'raw'/save name, column 2 =
%display name)
modelArray_disp = {'bcc-csm1-1','BCC-CSM1.1';...
                   'inmcm4','INM-CM4';...
                   'CESM1-0-4','CESM1.0.4';...
                   'CSIRO-Mk3-6-0','CSIRO Mk3.6.0';...
                   'IPSL','IPSL-CM5A-LR'};
             
%% Standard frequency bands to use
%Format for a data frequency: { [data frequency],{[n x 1 cell of long names
%for frequency bands]},[[n x 1 array of min/max *periods* for the frequency
%band]],{[n x 1 cell of 'IDs' for each frequency band - simple shorthands
%to use in function input. Should avoid spaces.]} }, for [n] frequency
%bands. Make sure the [data frequency] lines up with strcmp calls in
%Variability.m for new data frequencies added. 
freq_band_setup = {'day',{'< 5 days','5 - 30 days','30 - 365 days','> 1 year','all frequencies','2 - 15 days','15 - 90 days'},...
                         [0,5;5,30;30,365;365,Inf;0,Inf;2,15;15,90],...
                         {'HF','MF','LF','XF','Full','2to15','15to90'};...
                   'month',{'< 1 year','> 1 year','all frequencies','1 - 10 years','> 10 years'},...
                           [0,12;12,Inf;0,Inf;12,120;120,Inf],...
                           {'LF','XF','Full','XXF','XXXF'};...
                   'year',{'< 10 years',' > 10 years','all frequencies'},...
                          [0,10;10,Inf;0,Inf],...
                          {'XXF','XXXF','Full'};...
                   '3hour',{'< 5 days','5 - 30 days','30 - 365 days','> 1 year','all frequencies','2 - 15 days','15 - 90 days'},...
                           [0,40;40,240;240,365*8;365*8,Inf;0,Inf;16,15*8;15*8,720],...
                           {'HF','MF','LF','XF','Full','2to15','15to90'}};
                           
freq_band_setup_seasons = {'day',{'< 3 days','3 - 15 days','15 - 90 days'},...
							[0,3;3,15;15,90],...
							{'HF','MF','LF'};...
							'month',{'< 3 months'},...
							[0,3],...
							{'LF'}};

             
%% Frequency band name addendum
%For some graphing options, members of the following cell array will be
%followed by ' variability', to make it clearer ("< 1 year" by itself
%doesn't mean much....)
variability_ids = {'HF','MF','LF','XF','Full','XXF','XXXF','2to15','15to90'};

%% Save
var_list = who;
save('various_defaults','-v7.3',var_list{:})
disp(['various_defaults.mat saved, with variables: ',strjoin(var_list)])
