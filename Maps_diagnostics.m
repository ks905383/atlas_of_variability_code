function Maps_diagnostics(VarIndices,modelArray,varargin)
% MAPS_DIAGNOSTICS  Map various parameters, coll. known as 'pi diagnostics'
%
%   MAPS_DIAGNOSTICS(VarIndices,modelArray) creates a set of 4-panel maps
%   showing (by default) 1) PI mean values, 2) PI standard deviations, 3)
%   difference in mean values RCP8.5-PI, and 4) ratio of standard
%   deviations RCP8.5/PI for each of the variables defined by [VarIndices]
%   and models defined by [modelArray]. Maps are saved as .eps files in
%   [figure_dir]/[VarDomain]/PI_Diagnostics/.
%   
%   Warnings are thrown if graphs are not created for a specific
%   model/variable pair or if the calculations of the two different
%   experiments cover differing timeframes, but run continues. Warnings are
%   also thrown if a searched-for file does not exist (from NAME_CHARS
%   output, with all corresponding behavior).
%
%   Colorbar limits are defined by setting values across models using
%   MapLimits.m, taken from the resulting 'IntermodelStats' files. 
%   
%   Resulting figure is by default configured for 11" x 8" (792 x 576 pts)
%
%   Function requirements (on path): name_chars, var_chars, load_stddevs,
%   load_means, FileDomain, and (from MathWorks Community) structfind.
%
%   MAPS_DIAGNOSTICS(...,'[flag]',[params],...) modify program run as below:
%       'experiments',[cell]- manually set experiments to process. Default
%                             is {'piControl','rcp85'}; in general plots
%                             are 1) {1} mean values, 2) {1} std devs, 3)
%                             {2} - {1} mean values, 4) {2}/{1} std devs.
%       'season',[seas]     - show data separated by season, i.e.:
%           'JJA'               - June-July-August
%           'DJF'               - December-January-February
%       'close_graphs'      - close graphs in UI after saving
%       'pointsy',[int]     - change resultant image size (height, in
%                             points = 1/72th inch)
%       'pointsx',[int]     - change resultant image size (width, in 
%                             points = 1/72th inch)
%       'convert_png',[log] - set whether to save .png copy using gs
%                             (def: false)
%       'save_fig',[log]    - set whether to save .fig copy (def: false)
%       'replace',[log]     - set whether to replace existing files (def:
%                             true)
%       'color_split',[log],[int] - set whether to split colormap into
%                                   discrete colors (default: true, 32)
%       'testing'           - adds '_TEST' to filename
%       'save_log'          - exports log of attempted saves with
%                             success/failure status under
%                             [calc_dir]/MAPS_DIAGNOSTICS_log_[starttime].txt
%
%   Saving convention: 
%   [filevar]_[freq]_[model]_[exp2]_[exp1]_Diag_Maps(_[season]).eps   
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also NAME_CHARS, 
%   
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 11/26/2017

various_defaults = matfile('various_defaults.mat');

%Set clock (for logging purposes)
format shortG
startTimestamp = clock;

%Set display/saving settings
close_graphs = false;
pointsy = 576;
pointsx = 792;
save_log = false;
save_fig = false;
convert_png = false;
replace_files = true;
testing = false;
filename_add = [];
color_split = true;
colorbar_bands = 32;

%Set analysis settings
process_byseason = false; season = [];
folder_season = various_defaults.season_dir;
filename_add_save =[];

%Set experiments
expArray = {'piControl','rcp85'};

%Set behavior of optional function flags
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'experiments'}
                expArray = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'replace'}
                replace_files = varargin{in_idx+1};
            case {'close_graphs'}
                close_graphs = true;
            case {'pointsy'}
                pointsy = varargin{in_idx+1};
                if isa(varargin{in_idx+1},'numeric')==0 || varargin{in_idx+1} < 0
                    error('Invalid figure height points, must be a positive integer')
                end
            case {'pointsx'}
                pointsx = varargin{in_idx+1};
                if isa(varargin{in_idx+1},'numeric')==0 || varargin{in_idx+1} < 0
                    error('Invalid figure width points, must be a positive integer')
                end
            case {'save_log'}
                save_log = true;
            case {'convert_png'}
                convert_png = varargin{in_idx+1};
            case {'save_fig'}
                save_fig = varargin{in_idx+1};
            case {'testing'}
                testing = true;
                filename_add = '_TEST';
            case {'color_split'}
                color_split = varargin{in_idx+1};
                colorbar_bands = varargin{in_idx+2};
                if isa(varargin{in_idx+2},'numeric')==0 || varargin{in_idx+2} < 0
                    error('Invalid colorbar split, must be a positive integer')
                end
            case {'season'}
                process_byseason = true;
                season = varargin{in_idx+1};
        end
    end
end

num_models = length(modelArray);

%% 1 Setup Graph Characteristics
%Define colormap
ColorRtB = [0 0 0.3725;0 0 0.3804;0 0 0.3843;0 0 0.3922;0 0 0.3961;0 0 0.4039;0 0 0.4078;0 0 0.4157;0 0 0.4235;0 0 0.4275;0 0 0.4353;0 0 0.4392;0 0 0.4471;0 0 0.451;0 0 0.4588;0 0 0.4667;0 0 0.4706;0 0 0.4784;0 0 0.4824;0 0 0.4902;0 0 0.4941;0 0 0.502;0 0 0.5059;0 0 0.5137;0 0 0.5216;0 0 0.5255;0 0 0.5333;0 0 0.5373;0 0 0.5451;0 0 0.549;0 0 0.5569;0 0 0.5647;0 0 0.5686;0 0 0.5765;0 0 0.5804;0 0 0.5882;0 0 0.5922;0 0 0.6;0 0 0.6039;0 0 0.6118;0 0 0.6196;0 0 0.6235;0 0 0.6314;0 0 0.6353;0 0 0.6431;0 0 0.6471;0 0 0.6549;0 0 0.6627;0 0 0.6667;0 0 0.6745;0 0 0.6784;0 0 0.6863;0 0 0.6902;0 0 0.698;0 0 0.702;0 0 0.7098;0 0 0.7176;0 0 0.7216;0 0 0.7294;0 0 0.7333;0 0 0.7412;0 0 0.7451;0 0 0.7529;0 0 0.7608;0.0039 0.0039 0.7647;0.0196 0.0196 0.7686;0.0353 0.0353 0.7725;0.051 0.051 0.7765;0.0667 0.0667 0.7804;0.0824 0.0824 0.7843;0.098 0.098 0.7882;0.1137 0.1137 0.7922;0.1294 0.1294 0.7961;0.1451 0.1451 0.7961;0.1608 0.1608 0.8;0.1765 0.1765 0.8039;0.1922 0.1922 0.8078;0.2078 0.2078 0.8118;0.2235 0.2235 0.8157;0.2392 0.2392 0.8196;0.2549 0.2549 0.8235;0.2706 0.2706 0.8275;0.2863 0.2863 0.8314;0.302 0.302 0.8353;0.3176 0.3176 0.8392;0.3333 0.3333 0.8431;0.349 0.349 0.8471;0.3647 0.3647 0.851;0.3804 0.3804 0.8549;0.3961 0.3961 0.8588;0.4118 0.4118 0.8627;0.4275 0.4275 0.8627;0.4431 0.4431 0.8667;0.4588 0.4588 0.8706;0.4745 0.4745 0.8745;0.4902 0.4902 0.8784;0.5059 0.5059 0.8824;0.5216 0.5216 0.8863;0.5373 0.5373 0.8902;0.5529 0.5529 0.8941;0.5686 0.5686 0.898;0.5843 0.5843 0.902;0.6 0.6 0.9059;0.6157 0.6157 0.9098;0.6314 0.6314 0.9137;0.6471 0.6471 0.9176;0.6627 0.6627 0.9216;0.6784 0.6784 0.9255;0.6941 0.6941 0.9294;0.7098 0.7098 0.9333;0.7255 0.7255 0.9333;0.7412 0.7412 0.9373;0.7569 0.7569 0.9412;0.7725 0.7725 0.9451;0.7882 0.7882 0.949;0.8039 0.8039 0.9529;0.8196 0.8196 0.9569;0.8353 0.8353 0.9608;0.851 0.851 0.9647;0.8667 0.8667 0.9686;0.8824 0.8824 0.9725;0.898 0.898 0.9765;0.9137 0.9137 0.9804;0.9294 0.9294 0.9843;0.9451 0.9451 0.9882;0.9608 0.9608 0.9922;0.9765 0.9765 0.9961;0.9922 0.9922 1;1 0.9922 0.9922;0.9961 0.9765 0.9765;0.9922 0.9608 0.9608;0.9882 0.9451 0.9451;0.9843 0.9294 0.9294;0.9804 0.9137 0.9137;0.9765 0.898 0.898;0.9725 0.8824 0.8824;0.9686 0.8667 0.8667;0.9647 0.851 0.851;0.9608 0.8353 0.8353;0.9569 0.8196 0.8196;0.9529 0.8039 0.8039;0.949 0.7882 0.7882;0.9451 0.7725 0.7725;0.9412 0.7569 0.7569;0.9373 0.7412 0.7412;0.9333 0.7255 0.7255;0.9294 0.7098 0.7098;0.9255 0.6941 0.6941;0.9216 0.6784 0.6784;0.9176 0.6627 0.6627;0.9137 0.6471 0.6471;0.9098 0.6314 0.6314;0.9059 0.6157 0.6157;0.902 0.6 0.6;0.898 0.5843 0.5843;0.8941 0.5686 0.5686;0.8902 0.5529 0.5529;0.8863 0.5373 0.5373;0.8824 0.5216 0.5216;0.8784 0.5059 0.5059;0.8745 0.4902 0.4902;0.8745 0.4745 0.4745;0.8706 0.4588 0.4588;0.8667 0.4431 0.4431;0.8627 0.4275 0.4275;0.8588 0.4118 0.4118;0.8549 0.3961 0.3961;0.851 0.3804 0.3804;0.8471 0.3647 0.3647;0.8431 0.349 0.349;0.8392 0.3333 0.3333;0.8353 0.3176 0.3176;0.8314 0.302 0.302;0.8275 0.2863 0.2863;0.8235 0.2706 0.2706;0.8196 0.2549 0.2549;0.8157 0.2392 0.2392;0.8118 0.2235 0.2235;0.8078 0.2078 0.2078;0.8039 0.1922 0.1922;0.8 0.1765 0.1765;0.7961 0.1608 0.1608;0.7922 0.1451 0.1451;0.7882 0.1294 0.1294;0.7843 0.1137 0.1137;0.7804 0.098 0.098;0.7765 0.0824 0.0824;0.7725 0.0667 0.0667;0.7686 0.051 0.051;0.7647 0.0353 0.0353;0.7608 0.0196 0.0196;0.7569 0.0039 0.0039;0.7529 0 0;0.7451 0 0;0.7412 0 0;0.7333 0 0;0.7294 0 0;0.7216 0 0;0.7176 0 0;0.7098 0 0;0.7059 0 0;0.698 0 0;0.6941 0 0;0.6863 0 0;0.6824 0 0;0.6745 0 0;0.6667 0 0;0.6627 0 0;0.6549 0 0;0.651 0 0;0.6431 0 0;0.6392 0 0;0.6314 0 0;0.6275 0 0;0.6196 0 0;0.6157 0 0;0.6078 0 0;0.6039 0 0;0.5961 0 0;0.5882 0 0;0.5843 0 0;0.5765 0 0;0.5725 0 0;0.5647 0 0;0.5608 0 0;0.5529 0 0;0.549 0 0;0.5412 0 0;0.5373 0 0;0.5294 0 0;0.5255 0 0;0.5176 0 0;0.5098 0 0;0.5059 0 0;0.498 0 0;0.4941 0 0;0.4863 0 0;0.4824 0 0;0.4745 0 0;0.4706 0 0;0.4627 0 0;0.4588 0 0;0.451 0 0;0.4471 0 0;0.4392 0 0;0.4314 0 0;0.4275 0 0;0.4196 0 0;0.4157 0 0;0.4078 0 0;0.4039 0 0;0.3961 0 0;0.3922 0 0;0.3843 0 0;0.3804 0 0;0.3725 0 0];
%Set colormap
if color_split
    colormap_ratios = ColorRtB(1:(256/colorbar_bands):256,:);
    colormap_ratios = [colormap_ratios(2:colorbar_bands/2+1,:); colormap_ratios(colorbar_bands/2+1,:); colormap_ratios(colorbar_bands/2+2:end,:)];
    colormap_gray = flipud(gray(colorbar_bands));
else
    colormap_ratios = ColorRtB;
    colormap_gray = flipud(gray);
end
%Set Map Origin
map_origin = 155;
%Set projection
projection = 'bsam';

%Pre-allocate log array
if save_log
    complete_log = cell(1,length(VarIndices)*length(modelArray));
end

%Load directory / labeling defaults
various_defaults = matfile([various_defaults.code_dir,'various_defaults.mat']);
%Set experiment markers
exp_marker = cell(2,1);
for idx = 1:2
    try
        exp_marker(idx) = various_defaults.expArray_disp(find(cellfun(@(x) strcmp(expArray{idx},x),various_defaults.expArray_disp(:,1))),2); %#ok<FNDSB>
    catch
        exp_marker(idx) = expArray(idx);
    end
end

%% 2 Plot / Execute
for i = 1:length(VarIndices)
    %% 2.1 Get variable characteristics
    [~,filevarFN,freq,vardesc,units,~,clim_type,freqdesc,~] = var_chars(VarIndices(i));
    
    %% 2.2 Setup Data, Colors to Plot and Plot Labels
    %Load intermodel means/std devs for colorbar limit definitions
    IntermodelStats = cell(length(expArray),1); mean_idx = zeros(2,1); full_idx = zeros(2,1);
    for experiment = 1:2
        filenameIM = [various_defaults.proc_data_dir,'IntermodelStats/',filevarFN,freq,expArray{experiment},'_IntermodelStats.mat'];
        if exist(filenameIM,'file'); lim_file = matfile(filenameIM);
        else error('MapsDIAG:NoIntermodelStats',['The intermodel stats file for ',filevarFN,freq,expArray{experiment},...
                ' can not be found. Please check for its existence or generate it using the function MapLimits.m'])
        end
        IntermodelStats{experiment} = lim_file.IntermodelStats; clear lim_file
        %Find mean and full indices from IntermodelStats{1}
        mean_idx_tmp = structfind(IntermodelStats{experiment},'ID','mean');
        mean_idx(experiment) = mean_idx_tmp(1);
        full_idx(experiment) = structfind(IntermodelStats{experiment},'ID','std');
    end
    
    %Set colormaps and colormap limits for local mean map (different based
    %on variable characteristics)
    if strcmp(clim_type,'percent')==1
        colormap1=colormap_gray;
        if IntermodelStats{1}(mean_idx(1)).std/IntermodelStats(mean_idx(1)).mean+4*IntermodelStats{1}(mean_idx(1)).std<60;
            clim1=[0 IntermodelStats{1}(mean_idx(1)).std/IntermodelStats(mean_idx(1)).mean+4*IntermodelStats{1}(mean_idx(1)).std];
        else
            clim1=[0 100];
        end
    elseif strcmp(clim_type,'ratio')==1
        colormap1=colormap_gray;
        clim1=[0 IntermodelStats{1}(mean_idx(1)).mean+4*IntermodelStats{1}(mean_idx(1)).std];
    elseif strcmp(clim_type,'interval')==1
        colormap1=colormap_ratios;
        clim1=[-4*IntermodelStats{1}(mean_idx(1)).std 4*IntermodelStats{1}(mean_idx(1)).std];
    elseif strcmp(clim_type,'temp')==1
        colormap1=colormap_ratios;
        clim1=[273-3*IntermodelStats{1}(mean_idx(1)).std 273+3*IntermodelStats{1}(mean_idx(1)).std];
    end
    
    %Get positions for subplot axes
    position_sub = subplot_pos(4); %Base subplot positions
    pos_corx = [0.05 0 -0.01 0]; %Left column position correction 
    pos_corx2 = [0.04 0 -0.01 0]; %Right column position correction
    
    %Set up plot labels, what to plot
    PlotData = struct('Title',{['$\mu_{',exp_marker{1},'}$'],['$\sigma_{',exp_marker{1},'}$'],['$\mu_{',exp_marker{2},'}-\mu_{',exp_marker{1},'}$'],['$\log_{10}(\sigma_{',exp_marker{2},'}/\sigma_{',exp_marker{1},'})$']},...
                                    'colormap',{colormap1,colormap_gray,colormap_ratios,colormap_ratios},...
                                    'position',{position_sub(1,:)+pos_corx,position_sub(2,:)+pos_corx2,position_sub(3,:)+pos_corx,position_sub(4,:)+pos_corx2},...
                                    'colorbar_pos_adj',{[0,0.013,-0.4235,-0.035],[0.425,0.013,-0.4235,-0.035],[0,0.013,-0.4235,-0.035],[0.425,0.013,-0.4235,-0.035]},...
                                    'colorbar_label',{units,units,units,['\log_{10}(',units,' / ',units,')']},...
                                    'CLim',{clim1,[0 IntermodelStats{1}(full_idx(1)).mean+4*IntermodelStats{1}(full_idx(1)).std],[-abs(3*(IntermodelStats{2}(mean_idx(2)).mean - IntermodelStats{1}(mean_idx(1)).mean)) abs(3*(IntermodelStats{2}(mean_idx(2)).mean - IntermodelStats{1}(mean_idx(1)).mean))],[-0.5 0.5]});
    %clear full_idx, varname reused for different structs
    clear full_idx
                                
    %% 2.3 Plot for each model separately                            
    for j=1:num_models
        initial_vars = who;
        try
            %% 2.3.1 Model-Specific Graph Setup
            %Get filename characteristics
            if process_byseason
                dir_seas = [various_defaults.proc_data_dir,modelArray{j},folder_season];
                [~,run,strtyr,endyr] = name_chars(modelArray{j},expArray,filevarFN,freq,'directory',dir_seas,'file_end',['*Power_',season,'.mat']);
            else
                [~,run,strtyr,endyr] = name_chars(modelArray{j},expArray,filevarFN,freq);
            end
            
            %Set printed file name (to check existence for replacement)
            VarDomain = FileDomain(filevarFN);
            if process_byseason
                filename_add_save = ['_',season];
            end
            filenameS = [various_defaults.figure_dir,VarDomain,'/',filevarFN,freq,modelArray{j},'_',expArray{2},'_',expArray{1},'_',run{1},'Diag_Maps',filename_add_save];
            
            %% 2.3.2 Plot, if desired
            if replace_files || (~replace_files && ~exist([filenameS,'.eps'],'file'))
                %% 2.3.2.1 Load Data
                %Get file names/characteristics of files of inputted var/exp/model
                if ~process_byseason
                    [~,~,strtyr,~,~] = name_chars(modelArray{j},expArray,filevarFN,freq);
                    season_load = 'all';
                else
                    season_load = season;
                    [~,~,strtyr,~,~] = name_chars(modelArray{j},expArray,filevarFN,freq,'directory',dir_seas,'file_end',['*Power_',season,'.mat']);
                end
                
                %Load required data
                StdDevs(1).Raw = load_stddevs(VarIndices(i),modelArray{j},expArray{1},season_load,'start_year',strtyr{1},'season_folder',folder_season);
                Means(1).Raw = load_means(VarIndices(i),modelArray{j},expArray{1},season_load,'start_year',strtyr{1},'season_folder',folder_season);
                
                StdDevs(2).Raw = load_stddevs(VarIndices(i),modelArray{j},expArray{2},season_load,'start_year',strtyr{2},'season_folder',folder_season);
                Means(2).Raw = load_means(VarIndices(i),modelArray{j},expArray{2},season_load,'start_year',strtyr{2},'season_folder',folder_season);
                
                [~,lat,lon] = load_stddevs(VarIndices(i),modelArray{j},expArray{1},season_load,'start_year',strtyr{1},'season_folder',folder_season);
                
                %Find all frequencies data
                full_idx(1) = structfind(StdDevs(1).Raw,'ID','Full');
                full_idx(2) = structfind(StdDevs(2).Raw,'ID','Full');
                
                %% 2.3.2.2 Calculate Desired Quantities
                MeansDiff = Means(2).Raw - Means(1).Raw;
                StdDevRatio = StdDevs(2).Raw(full_idx(2)).Data./StdDevs(1).Raw(full_idx(1)).Data;
                
                if strcmp(modelArray{j},'CCSM4')==1 && (strcmp(filevarFN,'cllow')==1 || strcmp(filevarFN,'clmed')==1 || strcmp(filevarFN,'clhi')==1)
                    MeansDiff = 100*MeansDiff;
                    Means(1).Raw = 100*Means(1).Raw;
                    StdDevs(1).Raw(full_idx(1)).Data = 100*StdDevs(1).Raw(full_idx(1)).Data;
                end
                
                lStdDevRatio = log10(StdDevRatio);
                
                %Make sure lat/lon are double (otherwise mapping breaks)
                lat = double(lat);
                lon = double(lon);
                
                %Wrap data around 360 degrees so pcolorm doesn't produce dumb white line
                MeansDiff = [MeansDiff' MeansDiff(1,:)']';
                StdDevs(1).Raw(full_idx(1)).Data = [StdDevs(1).Raw(full_idx(1)).Data' StdDevs(1).Raw(full_idx(1)).Data(1,:)']';
                Means(1).Raw = [Means(1).Raw' Means(1).Raw(1,:)']';
                lStdDevRatio = [lStdDevRatio' lStdDevRatio(1,:)']';
                lon = [lon' 360]';
                lonWrapped = wrapTo360(lon);
                
                %Order data for plotting
                data_ordered = {Means(1).Raw,StdDevs(1).Raw(full_idx(1)).Data,MeansDiff,lStdDevRatio};
                for plot_idx = 1:4; PlotData(plot_idx).Data = data_ordered{plot_idx}; end; clear data_ordered;
                
                %% SET TO BETTER UNITS IF PRECIPITATION
                if ismember(VarIndices,[5,20,29])
                    if regexp(freq,'.[Dd]ay.')
                        unit_scale_constant = ((1/1000)*100*60*60*24);
                        units = '(cm / hr)';
                    else
                        unit_scale_constant = ((1/1000)*100*60*60*24*30);
                        units = '(cm / month)';
                    end
                    
                    clabel_tmp = {units,units,units,['\log_{10}(',units,' / ',units,')']};
                    [PlotData.('colorbar_label')]=clabel_tmp{:};
                    for plot_idx = 1:3
                        PlotData(plot_idx).Data = PlotData(plot_idx).Data*unit_scale_constant;
                        if j == 1
                            PlotData(plot_idx).CLim = PlotData(plot_idx).CLim*unit_scale_constant;
                        end
                    end
                end

                %% 2.3.2.3 Plot Data
                %Create Figures
                figure1 = figure;
                
                %Get World Map
                coast=load('coast');
                
                subplot_idx = cell(4,1);
                for plot_idx = 1:4
                    subplot_idx{plot_idx} = subplot('position',PlotData(plot_idx).position,'Visible','off','Parent',figure1,'CLim',PlotData(plot_idx).CLim,...
                        'FontSize',8,'XTick',[],'YTick',[],'DataAspectRatio',[1 1 1]);
                    ax = axesm(projection,'Origin',map_origin);
                    pcolorm(lat,lonWrapped,PlotData(plot_idx).Data.');
                    colormap(ax,PlotData(plot_idx).colormap);
                    c = colorbar(subplot_idx{plot_idx},'Position',PlotData(plot_idx).position+PlotData(plot_idx).colorbar_pos_adj);
                    c.Label.String = ['\fontsize{16}{0}$',PlotData(plot_idx).colorbar_label,'$']; c.Label.Interpreter = 'latex';
                    c.Label.Units = 'normalized';
                    geoshow(coast.lat,coast.long,'Color','black')
                    framem off; gridm off; axis off; mlabel off; plabel off;
                    title(['\fontsize{16}{0}',PlotData(plot_idx).Title],'Interpreter','latex')
                end
                
                %% 2.3.2.4 Annotate Figure
                %Set caption descriptions
                AnnotateString = [freqdesc,vardesc,char(10),modelArray{j},' ',season];
                
                %Set annotation
                annotation(gcf,'textbox',...
                    [0.05 0.01 0.9 0.0710172744721689],...
                    'String',{AnnotateString},...
                    'FitBoxToText','off',...
                    'FontSize',8,...
                    'EdgeColor','none');
                
                %Save as .fig, .eps, and .png
                figure1.PaperUnits = 'points';
                set(gcf,'PaperPositionMode','manual');
                figure1.PaperPosition = [0 0 pointsx pointsy];
                
                if testing
                    filenameS = [filenameS,filename_add]; %#ok<AGROW>
                end
                
                print(gcf,'-depsc2','-loose',filenameS);
                disp([filenameS,' saved'])
                if save_fig
                    saveas(figure1,filenameS)
                    disp([filenameS,'.fig saved']);
                end
                
                if convert_png
                    filenameSeps = [filenameS,'.eps'];
                    filenameSpng = [filenameS,'.png'];
                    system(strcat('gs -o -q -sDEVICE=png16m -dEPSCrop -r1200 -o', [filenameSpng,' ',filenameSeps]));
                    disp([filenameSpng,' saved']);
                end
                
                EndMsg=['Diagnostics Maps for ',modelArray{j},' ',filevarFN,freq,' Complete'];
                disp(EndMsg);
                
                if close_graphs
                    close(gcf);
                end
                
                %Store success message in log
                if save_log
                    save_log_msg = [filevarFN,freq,modelArray{j},' Created'];
                    complete_log{(i-1)*length(modelArray)+j} = save_log_msg;
                end
                
            else
                disp(['File for ',modelArray{j},' ',filevarFN,freq,num2str(strtyr{1}),'-',num2str(endyr{1}),' already exists, was not replaced']);
                if save_log
                    save_log_msg = [filevarFN,freq,modelArray{j},' already exists'];
                    complete_log{(i-1)*length(modelArray)+j} = save_log_msg;
                end
            end
            clearvars('-except',initial_vars{:});
        catch ME
            disp(ME)
            disp(ME.message)
            %Store warning in log
            if save_log
                save_log_wrn = ['****',filevarFN,freq,modelArray{j},' Not Complete****'];
                complete_log{(i-1)*length(modelArray)+j} = save_log_wrn;
            end
            warningMsg=['Diag Maps for ',modelArray{j},' ',filevarFN,freq,' Not Created!'];
            warning(warningMsg);
            clearvars('-except',initial_vars{:});
        end
    end
    EndMsg=['Diagnostics Maps for ',filevarFN,freq,' for All Models Complete'];
    disp(EndMsg);
    
end

%% 3 Export log as table, if desired
if save_log
    export_log = cell2table(complete_log','VariableNames',{'PIDiag1_attempted_executions'});
    writetable(export_log,[various_defaults.log_dir,'MAPS_DIAGNOSTICS_log_',num2str(startTimestamp(1)),'-',num2str(startTimestamp(2)),'-',num2str(startTimestamp(3)),'-',num2str(startTimestamp(4)),'-',num2str(startTimestamp(5)),'-',num2str(startTimestamp(6)),'.txt']);
end

end
