function panel_plots(VarIndices,modelArray,varargin)
% PANEL_PLOTS  plot log changes in variability vs log changes in 
%                       means
%
%   PANEL_PLOTS(VarIndices,modelArray) plots the log change in variability
%   for the LF, XF, and Full (usually corresponding to < 1 year, > 1 year,
%   and all frequencies, respectively) frequency bands in the data defined
%   by [VarIndices] and [modelArray] vs. log change in mean for each
%   variable in [VarIndices] and changes between RCP8.5 and pre-industrial
%   control experiment runs. The resulting figure has (num_freqs x
%   num_models) subplots, up to 6 models, beyond which new figures are
%   started until all models are plotted. The first principal component of
%   the data in each panel is calculated and graphed. An annotation in each
%   subplot gives Pearson's correlation coefficient, the slope of the first
%   principal component, and the percentage of points offscale (these are
%   also indicated on the map by points plotted on the corresponding axes
%   beyond which they are offscale). If more than one model is listed in
%   [modelArray], a set of cartoon axes below the panels shows the axes
%   labels; these are otherwise posted on the axes themselves. A caption is
%   printed, by default showing variable name and frequency, latitude bands
%   shown (0-90 absolute by default), figure number if more than one, and a
%   brief description of any modifications or clipping done through
%   optional flags, if relevant. The data is colorcoded by absolute
%   latitude, in a flipped jet-based colorscale with 18 colors, given by a
%   colorbar to the right of the panels. Graphs are saved as .eps files in
%   [figure_dir]/[VarDomain]. If only one model is selected, the file is
%   named after the model. Otherwise, the file contains the string
%   "ALLMODELS" in the filename model slot.
%
%   Data requirements: _sqrtPower, _LocalMeans, _StdDevCI (if
%   'remove_insig' flag used), for two experiments for the models listed in
%   [modelArray] for the variable given by [VarIndices]. If the
%   'scale_by_var' flag used, _LocalMeans is necessary for the variable
%   index following the flag as well. If the 'season' flag is used, the
%   _sqrtPower and _LocalMeans, etc., files must exist for the season
%   requested. If coloring points by a data or data ratio value and the
%   default intermodel coloring scheme is desired, a file _IntermodelStats
%   or _IntermodelRatioStats for the relevant variable and experiment(s) is
%   necessary (from function MapLimits). This code uses
%   [various_defaults.m] to set save directories and experiment display
%   names.
%
%   Function requirements (on path): clip_lat, clip_var_min, clip_var_max,
%   clip_ci, clip_region, data_colormap, load_means, load_stddevs, 
%   subplot_pos, FileDomain, name_chars, var_chars, and (from Mathworks 
%   community) dsn2fu and structfind.
%
%   PANEL_PLOTS(...,'[flag]',[params],...) modify program run 
%   as below:
%
%       Data management flags:
%       'experiments',[array]   - set scenarios to examine in a 2x1 cell
%                                 array, with comparison going from
%                                 [array]{2} to [array]{1} (def:
%                                 {'rcp85','piControl'}) (so, default
%                                 behavior is to show
%                                 log(RCP8.5/piControl) for data)
%       'clip_lat',[min],[max]  - restrict graph to specified (absolute)
%                                 latitude band (def: false, graphing 0-90)
%       'clip_region',[reg_id]  - restrict graph to region specified by
%                                 [reg_id] (either numerical index in 
%                                 Regions or string region abbreviation).
%                                 If no latitude clip chosen in addition, 
%                                 the filename does not show [0-90] as 
%                                 usual.
%       'clip_minval',[min],([log]) - clip points below a minimum mean
%                                     value [min] in second listed
%                                     experiment run in expArray (by
%                                     default PI Control). [log] determines
%                                     whether the thusly clipped points are
%                                     shown on the graph or not (they are
%                                     excluded from the PCA fit either
%                                     way). If no logical is provided, then
%                                     default behavior (not fitted points
%                                     are clipped/not displayed as well) is
%                                     kept. If true, the subplot annotation
%                                     additionally gives a percentage of
%                                     points displayed that are not fit
%                                     through. 
%       'clip_maxval',[max],([log]) - clip points below a maximum mean
%                                     value [max] in second listed
%                                     experiment run. Same behavior as
%                                     above.
%       'clip_insig'            - restrict graph to points whose 
%                                 variability ratio between given exps is 
%                                 'meaningfully' different from one as
%                                 determined in the bootstrap procedure
%                                 from VARIABILITY_CI
%       'hemisphere',[str]      - restrict graph to specified hemisphere
%                                 ('NH'/'SH')
%       'season',[str]          - restrict graph to specified season
%                                 ('DJF'/'JJA')
%       'scale_byvar',[var_idx] - scale data by local change in mean of a
%                                 variable given by index [var_idx] 
%       'show_freqs',[freqs]    - change frequencies shown. [freqs] can
%                                 be a cell array of either frequency band
%                                 sizes (i.e. {[0,5],[5,30],[30,365]}) or
%                                 IDs (i.e. {'HF','LF',...}) or a vector of
%                                 indices of the StdDev struct (discouraged
%                                 because the StdDev struct outputs from
%                                 Variability may not always be in a set
%                                 order)
%
%       Graph display flags: 
%       'close_graphs'          - close figures in the UI after saving
%       'visible',['off'/'on']  - set whether to open graphs in the UI at
%                                 all (def: 'on'; either option does not
%                                 affect other saving behavior of graphs or
%                                 stats)
%       'color_split',[int]     - set how many colors [int] the colormap is
%                                 divided into (default: 18, or 5 absolute
%                                 degrees latitude per color if default
%                                 'lat_abs' is chosen to color points by)
%       'pointsy',[int]         - change resultant image size (height, in
%                                 points = 1/72th inch)
%       'pointsx',[int]         - change resultant image size (width, in 
%                                 points = 1/72th inch)
%       'axes_cartoon',[log]    - set whether to show cartoon axes (def:
%                                 true for > 1 model, false for 1 model) 
%       'show_caption',[log]    - set whether to show caption with
%                                 variable and graph descriptors (def:
%                                 true)
%       'custom_caption',[str]  - set a custom caption
%       'show_title',[log]      - set whether to show title of figure (only
%                                 applies to single-model graphs; title in
%                                 that case is the model name + graphed
%                                 frequencies) (def: true)
%       'xlabel',[char]         - set xlabel (interpreter is latex)
%       'ylabel',[char]         - set ylabel (interpreter is latex)
%       'annotate_subplot',[log]- set whether to show annotation on subplot
%                                 showing slope of 1st principal component,
%                                 correlation coefficient, offscale points
%                                 percentage (if > 0), and displayed but
%                                 not fit through points (if > 0) (def:
%                                 true)
%       'show_legend',[log]     - set whether to show legend with fit (if
%                                 showing), data with lat subset, and large
%                                 dot at the mean value. Only relevant for
%                                 single-panel/single-subplot images (for
%                                 which default is true), does nothing for
%                                 other plotting schemas
%       'show_fit',[log]        - set whether to show the 1st principal
%                                 component on the graph (def: true)
%       'show_colorbar',[log]   - set whether to show colorbar (def: true,
%                                 but switches to false if colored by a
%                                 single color in color_by below)
%       'colorbar_label',[char],([char]) - manually set the colorbar label. 
%                                 If followed by a string 'tex', 'latex',
%                                 or 'none', then the text interpreter is
%                                 set as well. The default text interpreter
%                                 for the colorbar label (if none
%                                 specified) is 'latex' (MATLAB's limited 
%                                 use LaTeX ensemble)
%       'colorbar_range',[arr]  - manually set the colorbar range in a 2x1
%                                 numerical array (i.e. [-1 1]). By
%                                 default [0 90] (with 'lat_abs' as a
%                                 colorscheme), see below in 'color_by' - 
%                                 [color by type] description for default
%                                 color ranges for other coloring schemes.
%       'idx_count',[int]       - set whether to thin out graphed points
%                                 (by default, idx_count = 1 and all
%                                 unclipped points are plotted. For
%                                 idx_count ~= 1, only points
%                                 (1:idx_count:end) of the data arrays
%                                 after shuffling, clipping are plotted)
%       'point_size',[num]      - set size of scatterplot points (def: 4)
%       'marker_face_color',[char] - set whether to fill in the scatterplot
%                                    points ('flat') or not ('none',
%                                    default, showing rings as points to
%                                    allow for easier viewing of
%                                    overlapping points)
%
%       Graph coloring flags:
%       'color_by',____         - manually set what the colorbar should be
%                                 based on. By default, points are colored
%                                 by the absolute value of latitude. The
%                                 general form of these options is:
%       'color_by',[color by type],[exp idx(s)],[data type idx(s)],...
%       ([scale type]),([filename change])
%           
%           [color by type] = 'lat_abs': default, coloring each point by
%                                        the absolute value of its
%                                        latitude, in which case [exp
%                                        idx(s)] and [data type idx(s)] can
%                                        be omitted. Default colorbar range
%                                        [0 90].
%                           = 'lat':     coloring by actual (non absolute)
%                                        value of latitude, in which case
%                                        [exp idx(s)] and [data type
%                                        idx(s)] can be omitted. Default
%                                        colorbar range [-90 90].
%                           = 'data':    coloring each point by the value
%                                        of a different variable at the
%                                        same point. Default colorbar range
%                                        intermodel mean +/- 1 std, taken
%                                        from IntermodelStats (function
%                                        MapLimits).
%                           = 'data_ratio': coloring each point by the
%                                           value of the ratio of two
%                                           different variables at that
%                                           point, in which case [exp
%                                           idx(s)] and [data type idx(s)]
%                                           are both 2x1 cell arrays,
%                                           with coloring based on the
%                                           ratio of the 2nd element over
%                                           the 1st element. Automatically
%                                           switches to a Red to Gray to
%                                           Blue colormap. Default colorbar
%                                           range is intermodel mean +/- 1
%                                           std, taken from IntermodelStats
%                                           (function MapLimits).
%                           = 'color':   coloring all points the same
%                                        color, in which case [exp idx(s)]
%                                        is either a CSS3-supported color 
%                                        string (i.e. 'red') or an RGB
%                                        value. In this case, [data type
%                                        idx(s)] can be omitted. 
%           [exp idx(s)]        - can be either a string corresponding to
%                                 an experiment in expArray or a number
%                                 corresponding to an experiment's index in
%                                 the Data gathering structs below. If
%                                 'data_ratio', this is a 2x1 cell. If
%                                 'color', this is a color string or RGB
%                                 value. If 'lat', this is ignored.
%           [data type idx(s)]  - can be either a string corresponding to:
%                                   'mean': mean data
%                                [freq_id]: variability data, defined by 
%                                           freq_id just as in 'show_freqs' 
%                                 or an integer corresponding to a data
%                                 type's index in the Data gathering
%                                 structs below. If 'data_ratio', this is a
%                                 2x1 cell. If 'color' or 'lat', this is 
%                                 ignored.
%           ([scale type])      - optional, string corresponding to:
%                                   'raw': default, raw value of data or
%                                          data ratio
%                                   'log': natural log value of data or
%                                          data ratio
%                                 'log10': log base 10 value of data or
%                                          data ratio
%           ([filename change]) - optional, a filename add for this color
%                                 scheme. It MUST begin with '_' (i.e.
%                                 '_colorbypr'). Default filename adds are
%                                 as follows:
%                           for 'lat_abs': [] (none, default)
%                               for 'lat': '_colorbyactlat'
%                              for 'data': '_colorby[exp][freq_id]'
%                        for 'data_ratio': '_colorby[exps][freqs]ratio'
%                                           (if both exps or freqs are
%                                           equal, only one is added to
%                                           filename)
%                             for 'color': '_1color'
%       
%       Graph/Data saving flags:
%       'save_graphs',[log]     - set whether to save graphs (def: true)
%       'testing'               - add '_TEST' to filename
%       'save_png',[log]        - set whether to save .png copy
%                                 (def: false)
%       'save_stats',[log]      - set whether to save pca slope, corr
%                                 coefficients, graphed data, subsetting
%                                 parameters (from lat, region clips), and
%                                 frequency names in a file _stats by model
%                                 (def: false)
%       'filename_add',[str]    - add an extra custom suffix to the
%                                 filename
%       'save_dir',[str]        - change the directory to which figures are
%                                 outputted (by default, 'figure_dir' from
%                                 various_defaults.mat)
%
%   Saving convention: 
%   [filevar]_[freq]_'ALLMODELS'/[model]_[exp1]_[exp2]_Log_Scatter...
%                                              (_CI...)
%                                              (_Limited[min]-[max]...)
%                                              (_Local[varscale]Scaled...)
%                                              (_[region_id]...)
%                                               _[minlat]_[maxlat]...
%                                               _[shown frequency IDs]...
%                                              (_colorby[colorby desc.]...)
%                                               _[fig#]of[#]...
%                                               .eps
%
%   NOTE: this function is part of the Atlas of Variability code package,
%   but is in beta. 
%
%   NOTE: this code is designed to be relatively modular. To change loaded
%   data, see section "4.3.1.1 Load Data". To change how data is
%   manipulated, see "4.3.1.4 Calculate Quanitites to be Graphed and
%   Related Statistics". Correspondingly change axes labels under "1 Set
%   Display Settings" right after the introduction. This may have adverse
%   affects on color_by settings (they currently depend on the Data()
%   indices of mean and frequency-dependent data to stay constant - a
%   future version may fix this by tying 'mean' and '[ID]' calls in
%   color_by to track to freq_data and cons_data in section "4.3.1.2
%   Identify Frequency Bands to Plot").
%
%   NOTE: NaNs and Infs are removed from the dataset prior to calculations
%   of correlation and principal components.
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also  CLIP_LAT, CLIP_REGION, CLIP_VAR_MAX, CLIP_VAR_MIN, CLIP_CI,
%   LOAD_MEANS, LOAD_STDDEVS, VARIABILITY, VARIABILITY_CI, MAPLIMITS,
%   VARIOUS_DEFAULTS
%   
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 11/27/2017

%% 1 Set Display Settings
%Set axes labels
xlabel_string = '$\mathrm{log} (\mu_f/\mu_i)$';
ylabel_string = '$\mathrm{log} (\sigma_f/\sigma_i)$';

%Set display defaults (font sizes, position, etc.)
subp_pos_correction = [0.05 0.05 -0.05 -0.05];
DisplaySettings = struct('pointsx',{2560,2560,792,792},'pointsy',{1200,1200,576,576},...
                         'subp_position',{0,0,subplot_pos(3)+[subp_pos_correction; subp_pos_correction; subp_pos_correction; subp_pos_correction],subplot_pos(1)+[0.05 0.08 -0.05 -0.1]},...
                         'scaption_position',{[-0.95,0.95,0.1,0.1],[-0.95,0.95,0.1,0.1],[-0.65,0.7,0.1,0.1],[-0.65,0.7,0.1,0.1]},... %For axes caption(slope,corr,etc.) [0.24 0.81 0.04 0.04]
                         'caption_position',{[0.1 0.01 0.9 0.07],[0.1 0.01 0.9 0.07],[0.1 0.01 0.9 0.07],[0.1 0.01 0.9 0.07]},... %For figure caption ([var], [freq], etc.)
                         'axes_fontsize',{7,7,18,18},... %For x and y tickmarks and xlabels/ylabels (if not frequency or model) %CHANGED FROM 15 04/26
                         'title_fontsize',{11,11,18,18},... %For model and frequency labels
                         'caption_fontsize',{7,7,15,15},... %For axes caption (slope, corr, etc.)
                         'annotate_fontsize',{6,6,12,12},... %For figure caption ([var], [freq], etc.)
                         'colorbar_position',{[0.925 0.15 0.014 0.75],[0.925 0.15 0.014 0.75],[0.925 0.15 0.014 0.75],[0.85 0.15 0.025 0.75]},... %For colorbar (relative to colorbar axis)
                         'colorbar_fontsize',{13,13,18,18},... %For colorbar label %CHANGED FROM 15 04/26
                         'axes_cartoon',{true,true,false,false},... %Whether to print cartoon axes
                         'axes_cartoon_fontsize',{7,7,7,7},... %For cartoon axes
                         'line_width',{0.5,0.5,1.75,1.75},... %For all lines
                         'button_size',{12,12,30,30}); %For mean value button

%% 2 Set Defaults / Deal with Option Flags
%Gather directory / naming conventions
various_defaults = matfile('various_defaults.mat');

%Allow for non-cell inputs of a single model
if isa(modelArray,'char')
    modelArray = {modelArray};
end

%Clip/scaling defaults
expArray = {'rcp85','piControl'};
process_byhemisphere = false; hemisphere = [];
process_byseason = false; season = [];
folder_season = various_defaults.season_dir;
lat_clip = false; clip_lat_min = 0; clip_lat_max = 90;
clip_insig = false;
region_clip = false; region_id = [];
clip_minval = false; clip_minval_display = false;
clip_maxval = false; clip_maxval_display = false;
scale_byvar = false; scale_varidx = (11); 
weight_pca = false;

%Display defaults
axis_limits = [-1 1 -1 1];
show_caption = true; cust_caption = []; 
show_title = true;
show_freqs = {'LF','XF','Full'}; 
show_fit = true;
origin_cross = false;
origin_cross_dims = [-0.1 0.1 -0.1 0.1];
annotate_subplot = true;
show_legend = true;
close_graphs = false;
visible_set = 'on';
idx_count = 1;
point_size = 4;
marker_face_color = 'none';

%Color defaults
color_by = 'lat_abs';
color_data_scale = 'raw';
colorby_filename = [];
color_split = 18;
show_colorbar = true;
cust_clabel = [];
cinterpreter = 'latex';
cust_crange = [];

%Saving defaults
save_graphs = true;
save_stats = false;
save_png = false;
filename_add = [];
alt_save_dir = [];

%Set behavior of optional function flags
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            
            %Geographic/Clip/scaling function flags
            case {'experiments'}
                expArray = varargin{in_idx+1};
                varargin{in_idx+1} = 0; %(To allow for cell inputs, which would otherwise trip up switch)
            case {'hemisphere'}
                process_byhemisphere = true;
                hemisphere = varargin{in_idx+1};
            case {'season'}
                process_byseason = true;
                season = varargin{in_idx+1};
            case {'season_folder'}
                folder_season = varargin{in_idx+1};
            case {'clip_lat'}
                lat_clip = true;
                clip_lat_min = varargin{in_idx+1};
                clip_lat_max = varargin{in_idx+2};
                if ~isa(varargin{in_idx}+1,'numeric') || ~isa(varargin{in_idx}+2,'numeric')
                    error('Invalid latitude values; clip_lat should be followed by two integers, clip_lat_min and clip_lat_max');
                end
            case {'clip_insig'}
                clip_insig = true;
            case {'clip_region'}
                region_clip = true;
                region_id = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'clip_minval'}
                clip_minval = true;
                if length(varargin) >= in_idx+2
                    if islogical(varargin{in_idx+2})
                        clip_minval_display = varargin{in_idx+2};
                    end
                end
                %In case of precipitation, input in cm/yr (data in kg/m^2s)
                if VarIndices == 20 || VarIndices == 5 || VarIndices == 29
                    clip_minval_min = varargin{in_idx+1}/((1/1000)*100*60*60*24*365); %CHANGES UNITS INTO FROM cm/yr TO kg/m^2s (DATA UNITS)
                    clip_minval_disp = varargin{in_idx+1};
                else
                    clip_minval_min = varargin{in_idx+1}; clip_minval_disp = varargin{in_idx+1};
                end
            case {'clip_maxval'}
                clip_maxval = true;
                if length(varargin) >= in_idx+2
                    if islogical(varargin{in_idx+2})
                        clip_maxval_display = varargin{in_idx+2};
                    end
                end
                %In case of precipitation, input in cm/yr (data in kg/m^2s)
                if VarIndices == 20 || VarIndices == 5
                    clip_maxval_max = varargin{in_idx+1}/((1/1000)*100*60*60*24*365); %CHANGES UNITS INTO FROM cm/yr TO kg/m^2s (DATA UNITS)
                else
                    clip_maxval_max = varargin{in_idx+1};
                end
            case {'scale_byvar'}
                scale_byvar = true;
                scale_varidx = varargin{in_idx+1};
            case {'weight_pca'}
                weight_pca = varargin{in_idx+1};
            
            %Display Flags
            case {'axis_limits'}
                axis_limits = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'color_split'}
                color_split = varargin{in_idx+1};
                if isa(varargin{in_idx+1},'numeric')==0 || varargin{in_idx+1} < 0
                    error('Invalid colorbar split, must be a positive integer')
                end
            case {'pointsy'}
                [DisplaySettings(:).pointsy] = deal(varargin{in_idx+1}); %(in light of disp_idx not being id'ed till later, just set all .pointsy values to the input)
                if isa(varargin{in_idx+1},'numeric')==0 || varargin{in_idx+1} < 0
                    error('Invalid figure height points, must be a positive integer')
                end
            case {'pointsx'}
                [DisplaySettings(:).pointsx] = deal(varargin{in_idx+1});
                if isa(varargin{in_idx+1},'numeric')==0 || varargin{in_idx+1} < 0
                    error('Invalid figure width points, must be a positive integer')
                end
            case {'show_caption'}
                show_caption = varargin{in_idx+1};
            case {'custom_caption'}
                cust_caption = varargin{in_idx+1};
            case {'show_title'}
                show_title = varargin{in_idx+1};
            case {'xlabel'}
                xlabel_string = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'ylabel'}
                ylabel_string = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'axes_cartoon'}
                [DisplaySettings(:).axes_cartoon] = deal(varargin{in_idx+1});
            case {'show_freqs'}
                clear show_freqs
                show_freqs = varargin{in_idx+1};
                varargin{in_idx+1} = 0;  %(To allow for cell inputs, which would otherwise trip up switch)
            case {'show_fit'}
                show_fit = varargin{in_idx+1};
            case {'origin_cross'}
                if isa(varargin{in_idx+1},'logical')
                    origin_cross = varargin{in_idx+1};
                elseif isa(varargin{in_idx+1},'numeric')
                    origin_cross = true;
                    origin_cross_dims = varargin{in_idx+1};
                    varargin{in_idx+1} = 0;
                end
            case {'annotate_subplot'}
                annotate_subplot = varargin{in_idx+1};
            case {'show_legend'}
                show_legend = varargin{in_idx+1};
            case {'close_graphs'}
                close_graphs = true;
            case {'visible'}
                visible_set = varargin{in_idx+1};
            case {'show_colorbar'}
                show_colorbar = varargin{in_idx+1};
            case {'colorbar_label'}
                cust_clabel = varargin{in_idx+1};
                %if any(strcmp('tex',{'tex','latex','none'}))
                %   cinterpreter = varargin{in_idx+1}; 
                %end
            case {'color_range'}
                cust_crange = varargin{in_idx+1}; varargin{in_idx+1} = 0;
            case {'idx_count'} %Graphing (1:idx_count:end) of shuffled idxs
                idx_count = varargin{in_idx+1};
            case {'point_size'}
                point_size = varargin{in_idx+1};
            case {'marker_face_color'}
                marker_face_color = varargin{in_idx+1};
                
            %Color scheme flags
            case {'color_by'}
                color_by = varargin{in_idx+1};
                string_idx = [2,5]; %This counter eventually becomes the indices of all varargin inputs 2:5 from 'color_by' that are strings, for purposes of searching for filename, calc type below
                %Set what data to color by
                if strcmp(color_by,'data')
                    colorby_idx1 = varargin{in_idx+2}; varargin{in_idx+2} = 0; string_idx(1) = string_idx(1)+1;
                    colorby_idx2 = varargin{in_idx+3}; varargin{in_idx+3} = 0; string_idx(1) = string_idx(1)+1;
                    %Set experiment index in the Data struct (1st
                    %level) below (counter-intuitively, expArray{2}
                    %corresponds to Data(x,1) and vice-versa......)
                    if isa(colorby_idx1{1},'char') && any(~cellfun(@isempty,strfind(expArray,colorby_idx1{1})))
                        colorby_idx1_tmp(2) = find(cellfun(@isempty,strfind(expArray,colorby_idx1{1})));
                        %color_exp = [expArray{~cellfun(@isempty,strfind(expArray,colorby_idx1{1}))}];
                        color_exp = colorby_idx1{1};
                    end
                    %Set data type index in the Data struct (1st level)
                    if isa(colorby_idx2{1},'char') && strcmp(colorby_idx2{1},'mean')
                        colorby_idx1_tmp(1) = 2;
                    else
                        colorby_idx1_tmp(1) = 1;
                    end
                    %Allow identifying by frequency band ID, type by
                    %setting index in the Data struct (2nd level, .Raw 
                    %substruct)
                    if strcmp(colorby_idx2{1},'mean')
                           colorby_idx2_tmp(1) = 1; 
                           ctype = 'mean';
                    else colorby_idx2_tmp(1) = 0; ctype = 'variab'; %#ok<*SEPEX>
                    end
                    %(If 'std', and not 'mean, need actual data struct
                    %to determine index, is determined below)
                    
                    %Convert experiment, data type to linear index of the
                    %Data(). struct (1st level)
                    colorby_idx1 = sub2ind([2 2],colorby_idx1_tmp(1),colorby_idx1_tmp(2));
                    clear colorby_idx1_tmp;
                    %If both colorby data types are mean, convert data
                    %index to linear index of the Data().Raw(). struct (2nd
                    %level). If not, keep as a cell, need to identify
                    %frequency from data below. 
                    if ~ismember(0,colorby_idx2_tmp)%If mean (not freq-dependent)
                        colorby_idx2 = colorby_idx2_tmp; clear colorby_idx2_tmp; 
                        oper_desc{1} = '\mu'; %For default colorbar label
                    else
                        oper_desc{1} = '\sigma'; %For default colorbar label
                        %if both num/denom are freq-dependent, don't touch colorby_idx2, let the find frequency process figure it out below
                    end
                elseif strcmp(color_by,'data_ratio')
                    colorby_idx1 = varargin{in_idx+2}; varargin{in_idx+2} = 0; string_idx(1) = string_idx(1)+1;
                    colorby_idx2 = varargin{in_idx+3}; varargin{in_idx+3} = 0; string_idx(1) = string_idx(1)+1;
                    colorby_idx1_tmp = cell(2,1); colorby_idx2_tmp = zeros(2,1);
                    if isa(colorby_idx1{1},'char') 
                        %Set experiment index in the Data struct (1st
                        %level). Process iff both experiments are to be
                        %graphed / will have loaded data already (more
                        %general code not yet developed)
                        if find(~cellfun(@isempty,strfind(expArray,colorby_idx1{1}))) && find(~cellfun(@isempty,strfind(expArray,colorby_idx1{2})))
                            %Set experiment index in the Data struct (1st
                            %level) below (counter-intuitively, expArray{2}
                            %corresponds to Data(x,1) and vice-versa.....)
                            colorby_idx1_tmp{2}(1) = find(cellfun(@isempty,strfind(expArray,colorby_idx1{1}))); %So, as per commment above, look for the empty of the two
                            colorby_idx1_tmp{2}(2) = find(cellfun(@isempty,strfind(expArray,colorby_idx1{2})));
                            %Get names of the experiments of whom the
                            %ratios are taken for the eventual filename (if
                            %they are the same / the ratio is taken of two
                            %variables from the same experiment, only list
                            %one).
                            if find(cellfun(@isempty,strfind(expArray,colorby_idx1{1}))) ~= find(cellfun(@isempty,strfind(expArray,colorby_idx1{2})))
                                %Split them up into two, for colorbar
                                %labeling purposes later
                                color_exps{1} = expArray{~cellfun(@isempty,strfind(expArray,colorby_idx1{1}))};
                                color_exps{2} = expArray{~cellfun(@isempty,strfind(expArray,colorby_idx1{2}))};
                                color_exp = [color_exps{1},color_exps{2}];
                            else
                                color_exp = [expArray{~cellfun(@isempty,strfind(expArray,colorby_idx1{1}))}];
                                [color_exps{1:2}] = deal(expArray{~cellfun(@isempty,strfind(expArray,colorby_idx1{1}))});
                            end
                        else
                           error('LOG_PANEL:UnsupportedExp',['Not all color by experiments are ',...
                               'included in the standard data loading procedure, the code is not ',...
                               'yet set up to process coloring by non-loaded experimental data.'])
                        end
                        
                        %Set data type index in the Data struct (1st level)
                        if isa(colorby_idx2{1},'char') && strcmp(colorby_idx2{1},'mean')
                            colorby_idx1_tmp{1}(1) = 2; ctype{1} = 'mean';
                        else colorby_idx1_tmp{1}(1) = 1; ctype{1} = 'variab';
                        end
                        if isa(colorby_idx2{2},'char') && strcmp(colorby_idx2{2},'mean')
                            colorby_idx1_tmp{1}(2) = 2; ctype{2} = 'mean';
                        else colorby_idx1_tmp{1}(2) = 1; ctype{2} = 'variab';
                        end
                        
                        %Allow identifying by frequency band ID, type by
                        %setting index in the Data struct (2nd level, .Raw
                        %substruct). 1 == mean data (not
                        %frequency-dependent, which is what the struct
                        %indices usually capture), 0 identifies to be found
                        %frequency data. 
                        if strcmp(colorby_idx2{1},'mean')
                           colorby_idx2_tmp(1) = 1; 
                        else colorby_idx2_tmp(1) = 0; 
                        end
                        if strcmp(colorby_idx2{2},'mean')
                            colorby_idx2_tmp(2) = 1;
                        else colorby_idx2_tmp(2) = 0; 
                        end
                    end
                    %Convert experiment, data type to linear index of the
                    %Data(). struct (1st level)
                    colorby_idx1 = [sub2ind([2 2],colorby_idx1_tmp{1}(1),colorby_idx1_tmp{2}(1)) sub2ind([2 2],colorby_idx1_tmp{1}(2),colorby_idx1_tmp{2}(2))];
                    clear colorby_idx1_tmp;
                    %If both colorby data types are mean, convert data
                    %index to linear index of the Data().Raw(). struct (2nd
                    %level). If not, keep as a cell, need to identify
                    %frequency from data below. 
                    if ~ismember(0,colorby_idx2_tmp)%If both num/denom are mean (not freq-dependent)
                        colorby_idx2 = colorby_idx2_tmp; clear colorby_idx2_tmp; 
                        [oper_desc{1:2}] = deal('\mu'); %For default colorbar label
                    elseif any(colorby_idx2_tmp==1) %If one of num/denom is mean (not freq-dependent)
                        colorby_idx2{colorby_idx2_tmp==1} = 1;
                        %For default colorbar label, set description of
                        %operation
                        oper_desc{colorby_idx2_tmp==1} = '\mu'; %#ok<AGROW>
                        oper_desc{colorby_idx2_tmp==0} = '\sigma'; %#ok<AGROW>
                    else
                        [oper_desc{1:2}] = deal('\sigma'); %For default colorbar label
                        %if both num/denom are freq-dependent, don't touch colorby_idx2, let the find frequency process figure it out below
                    end
                    %For reference, colorby_idx1 is a 2x1 array, of which
                    %each element represents the linear index of Data() for
                    %each corresponding element. colorby_idx2 is the same,
                    %but with the linear index for Data.Raw() - but only if
                    %that index is already determined, otherwise it stays a
                    %cell.
                    
                elseif strcmp(color_by,'color')
                    colorby_idx1 = varargin{in_idx+2}; varargin{in_idx+2} = 0; 
                    %string_idx(1) = string_idx(1)+1;
                    string_idx = [2 2];
                    %Don't show colorbar if all one color
                    show_colorbar = false;
                elseif strcmp(color_by,'lat') || strcmp(color_by,'lat_abs')
                    string_idx = [2 2];
                end
                
                %Set data treatment for color scheme (raw or log)
                if find(~cellfun(@isempty,cellfun(@(x) strfind(x,'log10'),varargin(in_idx+string_idx(1):in_idx+string_idx(2)),'UniformOutput',0)))
                    color_data_scale = 'log10';
                elseif find(~cellfun(@isempty,cellfun(@(x) strfind(x,'log'),varargin(in_idx+string_idx(1):in_idx+string_idx(2)),'UniformOutput',0)))
                    color_data_scale = 'log';
                else %(by default)
                    color_data_scale = 'raw';
                end
                
                %Set def. filename change corresponding to color scheme as
                %able (for frequency-dependent, this will be updated later
                %to reflect frequency IDs from files)
                if isa(varargin{in_idx+string_idx(1)},'char') && all(cellfun(@(x) isa(x,'char'),varargin(in_idx+string_idx(1):in_idx+string_idx(2)))) && ~all(cellfun(@isempty,regexp(varargin(in_idx+string_idx(1):in_idx+string_idx(2)),'^_.')))
                    str_candidates = in_idx+string_idx(1):in_idx+string_idx(2);
                    colorby_filename = varargin{str_candidates(~cellfun(@isempty,regexp(varargin(in_idx+string_idx(1):in_idx+string_idx(2)),'_.')))};
                else
                   if strcmp(color_by,'lat_abs'); colorby_filename = [];end
                   if strcmp(color_by,'lat'); colorby_filename = '_colorbyactlat';end
                   if strcmp(color_by,'color'); colorby_filename = '_1color';end
                end
                
            %Saving Flags
            case {'save_graphs'}
                save_graphs = varargin{in_idx+1};
            case {'save_stats'}
                save_stats = true;
            case {'save_png'}
                save_png = varargin{in_idx+1};
            case {'testing'}
                filename_add = [filename_add,'_TEST']; %#ok<AGROW>
            case {'filename_add'}
                filename_add = [filename_add,varargin{in_idx+1}]; %#ok<AGROW>
            case {'save_dir'}
                alt_save_dir = varargin{in_idx+1};
        end
    end
end 

%% 3 Get Auxilliary Info
num_models = length(modelArray); % # models
num_freqs = length(show_freqs); % # frequency bands

%Set values for different formatting cases (indices of the struct with all 
%the display defaults DisplaySettings)
if length(modelArray) > 1
    if num_freqs>1; disp_idx = 1;
    else disp_idx = 2;
    end
    model_savename = 'ALLMODELS';
else
    if num_freqs>1; disp_idx = 3;
    else disp_idx = 4;
    end
    model_savename = modelArray{1};
end   

%Create subplot grid, based on number of models and frequency bands inputted
if num_models > 6 %(set up to do 6 models / figure, so creates more figures if more than 6 models)
    %Get number of figures, with 6 or fewer models per figure
    num_figures = floor(num_models/6);
    if mod(num_models,6)~=0
        num_figures = num_figures+1;
    end
    %Split up models by figure, in order of modelArray
    modelArray_sub = cell(1,num_figures);
    for figs = 1:num_figures-1
        modelArray_sub{figs} = modelArray(((figs-1)*6+1):figs*6);
    end
    clear figs
    modelArray_sub{num_figures} = modelArray((num_figures-1)*6+1:end);
    %Set subscripts to create 6 x num_frequencies graph grid
    model_vars_subscripts = (1:num_freqs*6);
    model_vars_subscripts = reshape(model_vars_subscripts,[6 num_freqs])';
    
else
    num_figures = 1;
    %Set subscripts for num_models x frequencies panel graph grid
    model_vars_subscripts = (1:num_freqs*num_models);
    model_vars_subscripts = reshape(model_vars_subscripts,[num_models num_freqs])';
    modelArray_sub{1} = modelArray;
end

%% 4 Plotting/Execution
for var_idx = VarIndices
    %% 4.1 Variable Setup
    %Define file variable identifiers
    [~,filevarFN,freq,vardesc,units,~,~,freqdesc,~] = var_chars(var_idx);
    if scale_byvar; [~,filevarFN_scale,~,vardesc_scale,~,~,~,freqdesc_scale,~] = var_chars(scale_varidx); end
    %Pre-allocated arrays for intermodel stats 
    y_data_mean = cell(length(show_freqs),length(modelArray)); x_data_mean = cell(length(show_freqs),length(modelArray));
    CorrCoeffs = cell(length(show_freqs),length(modelArray));
    slp = cell(length(show_freqs),length(modelArray)); mu = cell(length(show_freqs),length(modelArray));
    coeff = cell(length(show_freqs),length(modelArray)); explained = cell(length(show_freqs),length(modelArray));
    
    %% 4.2 Set base filename
    if ~isempty(alt_save_dir)
        filename_dir = alt_save_dir;
    else
        filename_dir = [various_defaults.figure_dir,FileDomain(filevarFN),'/'];
    end
    filename_base = [filevarFN,freq,model_savename,'_',expArray{1},'_',expArray{2},'_Log_Scatter'];
    if weight_pca; filename_base = strcat(filename_base,'_weightedpca');end
    if clip_insig; filename_base = strcat(filename_base,'_CI');end
    if clip_minval; filename_base = strcat(filename_base,'_Limited',num2str(clip_minval_disp));end
    if clip_minval && clip_maxval; filename_base = strcat(filename_base,'-',num2str(clip_maxval_disp)); end
    if clip_maxval && ~clip_minval; filename_base = strcat(filename_base,'_Limited',num2str(clip_maxval_disp));end
    if scale_byvar;filename_base = strcat(filename_base,'_LocalMean',upper(filevarFN_scale),'Scaled');end
    if region_clip; region_id_name = clip_region(region_id); 
        filename_base = strcat(filename_base,'_',region_id_name);end
    if lat_clip;filename_base = strcat(filename_base,'_',num2str(clip_lat_min),'_',num2str(clip_lat_max));
    elseif ~region_clip; filename_base = strcat(filename_base,'_',num2str(0),'_',num2str(90));
    end
    if isa(show_freqs{1},'char'); filename_base = strcat(filename_base,'_',show_freqs{:}); end
    if process_byhemisphere; filename_base = strcat(filename_base,'_',hemisphere); end
    if process_byseason; filename_base = strcat(filename_base,'_',season); end
    
    %% 4.3 Create Plots
    for f = 1:num_figures
        %% 4.3.1 Set up models for current figure
        modelArray = modelArray_sub{f};
        num_models = length(modelArray);
        
        %Create figure with a print size determined by pointsx, pointsy
        figure1 = figure('visible',visible_set,'PaperUnits','points','PaperPositionMode','manual','PaperPosition',[0 0 DisplaySettings(disp_idx).pointsx DisplaySettings(disp_idx).pointsy]);
        
        %% 4.3.2 Plot by Model
        for j=1:num_models
            %% 4.3.2.1 Load Data
            %Get file names/characteristics of files of inputted var/exp/model
            if ~process_byseason
                season_load = 'all';
                [~,~,strtyr,~,~] = name_chars(modelArray{j},expArray,filevarFN,freq);
            else
                season_load = season;
                folder = folder_season; dir_seas = [various_defaults.proc_data_dir,modelArray{j},folder];
                [~,~,strtyr,~,~] = name_chars(modelArray{j},expArray,filevarFN,freq,'directory',dir_seas,'file_end',['*Power_',season,'.mat']);
            end
            
            %Load required data for first variable
            Data(1,1).Raw = load_stddevs(var_idx,modelArray{j},expArray{2},season_load,'start_year',strtyr{2},'season_folder',folder_season);
            Data(2,1).Raw = load_stddevs(var_idx,modelArray{j},expArray{2},season_load,'start_year',strtyr{2},'season_folder',folder_season);
            Data(2,1).Raw = Data(2,1).Raw(structfind(Data(2,1).Raw,'ID','mean')).Data;
            
            %Load required data for second variable
            Data(1,2).Raw = load_stddevs(var_idx,modelArray{j},expArray{1},season_load,'start_year',strtyr{1},'season_folder',folder_season);
            Data(2,2).Raw = load_stddevs(var_idx,modelArray{j},expArray{1},season_load,'start_year',strtyr{1},'season_folder',folder_season);
            Data(2,2).Raw = Data(2,2).Raw(structfind(Data(2,2).Raw,'ID','mean')).Data;
            
            %Load area scaling weights for mean calculations
            filename_w = strcat(various_defaults.code_dir,'GridWeights/',modelArray{j},'_GridWeights');
            load(filename_w); GridWeights_tmp = GridWeights; clear GridWeights
            
            %Load Mean of a different variable, if scaling by var Change
            if scale_byvar
                LocalMeans_scale{1} = load_means(scale_varidx,modelArray{j},expArray{2},'start_year',strtyr{2});
                LocalMeans_scale{2} = load_means(scale_varidx,modelArray{j},expArray{1},'start_year',strtyr{1});
                
                local_change_scale_tmp = LocalMeans_scale{2}-LocalMeans_scale{1};
            end
        
            %Get geo array
            [lat,lon] = load_stddevs(var_idx,modelArray{j},expArray{1},season_load,'start_year',strtyr{1},'season_folder',folder_season);
            
            %% 4.3.2.2 Identify Frequency Bands to Plot
            struct_idxs = zeros(numel(Data),1);
            for data_idx = 1:numel(Data) %Find which loaded data is stored in structs (and therefore is frequency band-dependent)
                if isa(Data(data_idx).Raw,'struct'); struct_idxs(data_idx) = 1; end
            end
            freq_data = find(struct_idxs);
            cons_data = find(~struct_idxs);
            
            %If any data is frequency band-dependent, determine the indices
            %of the desired bands in each struct
            if any(struct_idxs)
                idx = cell(length(freq_data),length(show_freqs));
                %Determine frequency indices for data to be graphed
                for freq_data_idx = 1:length(freq_data)
                    for band = 1:length(show_freqs)
                        if isa(show_freqs,'numeric') %Determine frequency bands by just indexes (discouraged)
                            if max(show_freqs) <= length(Data(freq_data(freq_data_idx)).Raw) %Make sure the desired number of frequency bands exists in the data
                                idx(freq_data_idx,:) = num2cell(show_freqs);
                            else
                                idx(freq_data_idx,1:length(Data(freq_data(freq_data_idx)).Raw)) = num2cell(show_freqs(1:length(Data(freq_data(freq_data_idx)))));
                                idx(freq_data_idx,length(Data(freq_data(freq_data_idx)).Raw)+1:end) = [];
                            end
                        elseif isa(show_freqs,'cell') %Determine frequency bands by size or by ID
                            if isa(show_freqs{band},'numeric') %Determined by size
                                idx{freq_data_idx,band} = structfind(Data(freq_data(freq_data_idx)).Raw,'Size',show_freqs{band});
                            elseif isa(show_freqs{band},'char') %Determined by ID
                                idx{freq_data_idx,band} = structfind(Data(freq_data(freq_data_idx)).Raw,'ID',show_freqs{band});
                            end
                        end
                    end
                    
                    %Determine frequency indices for data to be colored by
                    %(if not already determined/converted to matrix above)
                    if (strcmp(color_by,'data') || strcmp(color_by,'data_ratio')) && isa(colorby_idx2,'cell')
                        for color_el = find(cellfun(@(x) isa(x,'char'),colorby_idx2)) %For number of frequency-dependent color components (non-freq dependent == {1}, as determined in flag/option section)
                            if isa(show_freqs{color_el},'numeric') %Determined by size
                                colorby_idx2{color_el} = structfind(Data(freq_data(freq_data_idx)).Raw,'Size',colorby_idx2{color_el});
                            elseif isa(show_freqs{color_el},'char') %Determined by ID
                                colorby_idx2{color_el} = structfind(Data(freq_data(freq_data_idx)).Raw,'ID',colorby_idx2{color_el});
                            end
                        end
                        %Convert from cell to linear indices of
                        %Data().Raw(). (freq band data, 2nd level) 
                        colorby_idx2 = cell2mat(colorby_idx2);
                    end
                end
                if any(isempty(idx{freq_data_idx,band}))
                    error('Plots_Log_Panel:NoFreqAvail',['The frequency bands ',show_freqs{:},' are not all available in both scenario runs for ',modelArray{j},' ',filevarFN])
                    %INSERT WARNING / COPING MECHANISM HERE
                    %freq_ids = something something intersect
                    %Based on which columns are not empty
                else
                    freq_ids = cell2mat(idx); clear idx
                end
            else
                %freq_ids = ones(1,2);
                freq_ids = [];
            end
            
            %Get final number of frequencies to plot
            num_freqs = size(freq_ids,2);
            
            %Isolate colorby_data. If colorby_idx2 == 1, then the data is
            %mean, and not in the same struct format as the variability
            %data, and has to be dealt with slightly differently (.Raw is
            %the last level in this struct in this case)
            if strcmp(color_by,'data') || strcmp(color_by,'data_ratio')
                if colorby_idx2(1) ~= 1; colorby_data{1} = Data(colorby_idx1(1)).Raw(colorby_idx2(1)).Data;
                else colorby_data{1} = Data(colorby_idx1(1)).Raw; end
                if strcmp(color_by,'data_ratio')
                    if colorby_idx2(2)~=1; colorby_data{2} = Data(colorby_idx1(2)).Raw(colorby_idx2(2)).Data;
                    else colorby_data{2} = Data(colorby_idx1(2)).Raw;
                    end
                end
            end
            
            %% 4.3.2.3 Subset Data
            if any(struct_idxs) %If frequency-dependent data, allow for possibilitiy of frequency-dependent clipping
                idx_display_clip = cell(num_freqs,1); 
                idx_fit_clip = cell(num_freqs,1); %If certain points are to be plotted, but not fit through
                area_points = cell(num_freqs,1); %Number of base points in geographic region
            else
                idx_display_clip = {[]};
                idx_fit_clip = {[]};
                area_points = {[]};
            end
            
            if process_byhemisphere; [idx_display_clip,~,~,~,~] = clip_lat(0,90,lat,lon,idx_display_clip,'hemisphere',hemisphere);
            end

            if lat_clip; [idx_display_clip,~,~,~,~] = clip_lat(clip_lat_min,clip_lat_max,lat,lon,idx_display_clip);
            end
            
            if region_clip; [idx_display_clip,~,~,~,~,region_id] = clip_region(region_id,modelArray{j},lat,lon,idx_display_clip);
            end
            
            %Get number of points in the geographic area isolated
            for vars = 1:num_freqs; area_points{vars} = length(lon)*length(lat)-length(idx_display_clip{vars}); end
            
            if clip_insig; [idx_display_clip,~,~,num_ciclipped] = clip_ci(var_idx,modelArray{j},expArray,idx_display_clip,'use_freqs',{Data(freq_data(1)).Raw(freq_ids(1,:)).ID});
            end
            
            %After all clip behaviors that solely affect idx_display_clip,
            %idx_fit_clip = idx_display_clip U [idxs of points not in fit, but shown]
            if clip_minval
                if ~clip_minval_display
                    [idx_display_clip,~,~,~] = clip_var_min(clip_minval_min,load_means(var_idx,modelArray{j},expArray{2},season_load,'start_year',strtyr{2},'season_folder',folder_season),idx_display_clip);
                else [~,idx_fit_clip,~,~] = clip_var_min(clip_minval_min,load_means(var_idx,modelArray{j},expArray{2},season_load,'start_year',strtyr{2},'season_folder',folder_season),idx_display_clip);
                end
            end
            
            if clip_maxval
                if ~clip_maxval_display
                    [idx_display_clip,~,~,~] = clip_var_max(clip_maxval_max,load_means(var_idx,modelArray{j},expArray{2},season_load,'start_year',strtyr{2},'season_folder',folder_season),idx_display_clip);
                else [~,idx_fit_clip,~,~] = clip_var_max(clip_maxval_max,load_means(var_idx,modelArray{j},expArray{2},season_load,'start_year',strtyr{2},'season_folder',folder_season),idx_display_clip);
                end
            end
            
            %Throw error if too much is clipped
            if any(cellfun(@length,idx_display_clip)==length(lat)*length(lon))
                error('PANEL_PLOTS:TooMuchClipped','no display of points possible; all points are clipped, please choose less restrictive clipping limits');
            end
            
            %Clip points determined through the process above. 
            %freq_ids = [# frequency-dependent variables, in order of
            %struct_idxs] x num_freqs matrix giving index of each frequency
            %within the Data.Raw() struct for each Data() that is
            %frequency-dependent for the data of that frequency band
            %freq_data is find(struct_idxs) (based on the thought that all
            %struct variables are frequency-dependent, by file convention)
            if any(struct_idxs) %If any data is frequency-dependent
               cons_data_tmp = cell(length(cons_data),length(Data(freq_data(1)).Raw));
               for vars = 1:num_freqs
                   for k = 1:length(freq_data)
                       [~,k2] = ind2sub(size(Data),freq_data(k)); %Get index of relevant frequency band
                       Data(freq_data(k)).Raw(freq_ids(k2,vars)).Data(idx_display_clip{vars}) = []; %Remove clipped indices determined above
                   end
                   for k = 1:length(cons_data) %Constant (not frequency-dependent) data is replicated num_freqs (length of freq_ids) times, since subsetting above may depend on frequency band-dependent values
                       [~,k2] = ind2sub(size(Data),cons_data(k)); %Get index of relevant constant data
                       cons_data_tmp{k,freq_ids(k2,vars)} = Data(cons_data(k)).Raw; %Deal out constant data num_freqs times
                       cons_data_tmp{k,freq_ids(k2,vars)}(idx_display_clip{vars}) = []; %Remove clipped indices determined above
                   end
               end
               for k = 1:length(cons_data)
                   %Turn the constant data into a struct to share data
                   %structure of above/store data in a standardized fashion
                   Data(cons_data(k)).Raw = struct('Data',cons_data_tmp(k,:)); 
                   Data(cons_data(k)).Raw(1).ID = 'mean'; %For colorscaling description's sake
               end
               if scale_byvar %Same treatment as constant data above, for scale variable, if desired
                   [local_change_scale{1:num_freqs}] = deal(local_change_scale_tmp); clear local_change_scale_tmp
                   for vars = 1:num_freqs
                       local_change_scale{vars}(idx_display_clip{vars}) = [];
                   end
               end
               clear cons_data_tmp;
            else %If no data is frequency-dependent, process the same for all components of Data
                for k = 1:numel(Data) %To make sure data is stored in a standardized fashion for selection below
                    cons_data_tmp = Data(k).Raw(idx_display_clip{1});
                    Data(k).Raw.ID = 'mean'; %For colorscaling description's sake
                    Data(k).Raw.Data = cons_data_tmp;
                    Data(k).Raw.Data(idx_display_clip{1}) = [];
                    if scale_byvar
                        local_change_scale(idx_display_clip{1}) = [];
                    end
                end
                clear cons_data_tmp
            end
            
            %Clip points on grid weights for mean calculations
            [GridWeights{1:num_freqs}] = deal(GridWeights_tmp); clear GridWeights_tmp
            for vars = 1:num_freqs; GridWeights{vars}(idx_display_clip{vars}) = []; end
            
            %Set points to be fit through (potentially different from
            %points displayed)
            if any(cellfun(@isempty,idx_fit_clip)) %If no special subsetting for points to be fit through, set equal to points displayed
                idx_fit_clip = idx_display_clip;
            end
            
            %If number of points clipped is the same length as number of
            %pixels (nothing to be plotted), throw error
            if any(cellfun(@length,idx_display_clip)==length(lat)*length(lon))
                error('PANEL_PLOTS:TooMuchClipped','no fit calculations possible; all points are clipped, please choose less restrictive clipping limits');
            end
            
            %% 4.3.2.4 Calculate Quantities to be Graphed and Related Statistics
            y_data = cell(num_freqs,1); x_data = cell(num_freqs,1); shuffled_idxs = cell(num_freqs,1); area_weights = cell(num_freqs,1);
            %Get a unified random permutation of indices to randomize
            %plotting order (in order to avoid systemic hiding of trends
            %underneath other points)
            if size(lon,2) == 1
                shuffled_idxs_tot = randperm(length(lon)*length(lat) - max(cellfun(@length,idx_display_clip)));
            else
                shuffled_idxs_tot = randperm(size(lon,1)*size(lon,2) - max(cellfun(@length,idx_display_clip)));
            end
            
            for vars = 1:num_freqs
                %Get frequency-band specific shuffled idxs (in case of
                %different clip behavior for different bands, leading to a 
                %different number of pixels to graph between bands)
                shuffled_idxs{vars} = shuffled_idxs_tot(shuffled_idxs_tot <= length(lon)*length(lat)-length(idx_display_clip{vars}));
                %Set subset of points, to simplify/save space, if desired
                shuffled_idxs{vars} = shuffled_idxs{vars}(1:idx_count:end);
                
                %Set data to be graphed
                y_data{vars} = log(Data(1,2).Raw(freq_ids(2,vars)).Data(shuffled_idxs{vars})./Data(1,1).Raw(freq_ids(1,vars)).Data(shuffled_idxs{vars}));
                x_data{vars} = log(Data(2,2).Raw(freq_ids(2,vars)).Data(shuffled_idxs{vars})./Data(2,1).Raw(freq_ids(1,vars)).Data(shuffled_idxs{vars}));
                
                %Scale data if desired
                if scale_byvar
                    y_data{vars} = y_data{vars}./local_change_scale{vars}(shuffled_idxs{vars});
                    x_data{vars} = x_data{vars}./local_change_scale{vars}(shuffled_idxs{vars});
                end
                
                %Remove Infs and NaNs (can arise through IEEE 0/0, 1/0
                %standards, Infs throw off PCA calculations below)
                rm_idxs = unique([find(isinf(y_data{vars})),find(isinf(x_data{vars})),find(isnan(y_data{vars})),find(isnan(x_data{vars}))]);
                y_data{vars}(rm_idxs) = []; x_data{vars}(rm_idxs) = [];
                shuffled_idxs{vars}(rm_idxs) = [];
                
                %Calculate correlation between graphed data and each axis'
                %(area-weighted) mean
                CorrCoeffs{vars,(f-1)*6+j} = corrcoef(x_data{vars},y_data{vars},'rows','complete');
                y_data_mean{vars,(f-1)*6+j} = (nansum(y_data{vars}.*(GridWeights{vars}(shuffled_idxs{vars}))))./nansum(GridWeights{vars}(shuffled_idxs{vars}));
                x_data_mean{vars,(f-1)*6+j} = (nansum(x_data{vars}.*(GridWeights{vars}(shuffled_idxs{vars}))))./nansum(GridWeights{vars}(shuffled_idxs{vars}));
                
                %Set area weights for pca, if desired, using GridWeights
                if weight_pca;area_weights{vars} = GridWeights{vars}(shuffled_idxs{vars}); else area_weights{vars} = ones(length(shuffled_idxs{vars}),1)/length(shuffled_idxs{vars}); end
            end
            
            %% 4.3.2.5 Set Colors of Points
            
            %Get lon/lat grid for colormapping
            if size(lon,2) == 1
                [~, latyy] = ndgrid(lon, lat);
            else
                latyy = lat;
            end
            %Set data to color by and associated colormap ranges
            if strcmp(color_by,'lat_abs')
                color_data = abs(latyy);
                %Set default colormap range
                def_crange = [0 90];
                %Set defaul colorbar label
                def_clabel = '\mathrm{Latitude} \ \mathrm{(absolute)}';
            elseif strcmp(color_by,'lat')
                color_data = latyy;
                %Set default colormap range
                def_crange = [-90 90];
                %Set defaul colorbar label
                def_clabel = '\mathrm{Latitude}';
            elseif strcmp(color_by,'data')
                %color_data = Data(colorby_idx1).Raw(colorby_idx2).Data;
                color_data = colorby_data{1};
                color_freq = Data(colorby_idx1).Raw(colorby_idx2).ID;
                if isempty(colorby_filename); colorby_filename = ['_colorby',color_exp,color_freq];
                end
                %Set default colormap range
                def_crange = 'def';
                %Make the subscript of the colorbar label easier to read
                if strcmp(color_exp,'rcp85'); color_exp = 'RCP'; 
                elseif strcmp(color_exp,'piControl'); color_exp = 'PI';
                elseif strcmp(color_exp,'historical'); color_exp = 'HIST';
                end
                %Set defaul colorbar label
                def_clabel = [oper_desc{1},'(',color_freq,')_{',color_exp,'}'];
            elseif strcmp(color_by,'data_ratio') %MAKE SURE ORDER OF INDICES IS CORRECT
                %color_data = Data(colorby_idx1(1)).Raw(colorby_idx2(1)).Data./Data(colorby_idx1(2)).Raw(colorby_idx2(2)).Data;
                color_data = colorby_data{1}./colorby_data{2};
                if colorby_idx2(1) == colorby_idx2(2)
                    color_freq = Data(colorby_idx1(1)).Raw(colorby_idx2(1)).ID;
                else
                    color_freq = [Data(colorby_idx1(1)).Raw(colorby_idx2(1)).ID,Data(colorby_idx1(2)).Raw(colorby_idx2(2)).ID];
                end
                if isempty(colorby_filename); colorby_filename = ['_colorby',color_exp,color_freq,'ratio'];
                end
                %Set default colormap range
                def_crange = 'def';
                %Make the subscript of the colorbar label easier to read
                for ridx = 1:2
                    if strcmp(color_exps{ridx},'rcp85'); color_exps{ridx} = 'RCP'; 
                    elseif strcmp(color_exps{ridx},'piControl'); color_exps{ridx} = 'PI';
                    elseif strcmp(color_exps{ridx},'historical'); color_exps{ridx} = 'HIST';
                    end
                end
                %Set default colorbar label
                def_clabel = [oper_desc{1},'(',Data(colorby_idx1(1)).Raw(colorby_idx2(1)).ID,')_{',color_exps{1},'} / ',...
                    oper_desc{2},'(',Data(colorby_idx1(2)).Raw(colorby_idx2(2)).ID,')_{',color_exps{2},'}'];
            elseif strcmp(color_by,'color')
                color_data = ones(length(lon)*length(lat),1);
                %Set default colormap range
                def_crange = [0 0];
                %Set defaul colorbar label
                def_clabel = '';
                %Supress data_colormap warning on too few colors, since it
                %warns about something that is done by design
                warning('off','DATA_COLORMAP:TooFewColors')
            end
            
            %Get colorbar range, either default, calculated from
            %IntermodelStats, or manual (set by flag)
            color_warn = false; %Set warning to off by default, is otherwise triggered under certain conditions below
            if ~isempty(cust_crange)
                c_range = cust_crange;
            else %If no custom color range used
                %Attempt to use intermodel stats to set an intermodel
                %applicable color range
                if strcmp(color_by,'data_ratio')
                    if exist([various_defaults.proc_data_dir,'IntermodelStats/',filevarFN,freq,expArray{colorby_idx1(1)},'_',expArray{colorby_idx1(2)},'_IntermodelRatioStats.mat'],'file')
                        load([various_defaults.proc_data_dir,'IntermodelStats/',filevarFN,freq,expArray{colorby_idx1(1)},'_',expArray{colorby_idx1(2)},'_IntermodelRatioStats.mat'])
                        %Find relevant statistic in Intermodel Stats
                        if strcmp(ctype,'mean');cfreq_idx = structfind(Intermodel_Ratio_Stats,'ID','mean');
                        else cfreq_idx = structfind(Intermodel_Ratio_Stats,'ID',Data(colorby_idx1(1)).Raw(colorby_idx1(1)).ID);
                        end
                        %Set c_range to mean +/- 1 standard deviation
                        c_range = [Intermodel_Ratio_Stats(cfreq_idx).mean-Intermodel_Ratio_Stats(cfreq_idx).std Intermodel_Ratio_Stats(cfreq_idx).mean+Intermodel_Ratio_Stats(cfreq_idx).std];
                        %Adjust to log/log10, if color data is log/log10
                        if strcmp(color_data_scale,'log'); c_range = log(c_range);
                        elseif strcmp(color_data_scale,'log10'); c_range = log10(c_range);
                        end
                    else
                        color_warn = true;
                        warning('LOG_PANEL:NoIntermStats',['No intermodel stats found for ',filevarFN,freq,...
                            ' ',expArray{colorby_idx(1)},' and ',expArray{colorby_idx1(2)},...
                            ', reverting to default color range (from min to max of ',...
                            'color by dataset, determined individually by model'])
                        c_range = def_crange;
                    end
                elseif strcmp(color_by,'data')
                    if exist([various_defaults.proc_data_dir,'IntermodelStats/',filevarFN,freq,expArray{colorby_idx1(1)},'_IntermodelStats.mat'],'file')
                        %Load intermodel stats
                        load([various_defaults.proc_data_dir,'IntermodelStats/',filevarFN,freq,expArray{colorby_idx1(1)},'_IntermodelStats.mat'])
                        %Find relevant statistic in Intermodel Stats
                        if strcmp(ctype,'mean');cfreq_idx = structfind(IntermodelStats,'ID','mean');
                        else cfreq_idx = structfind(IntermodelStats,'ID',Data(colorby_idx1(1)).Raw(colorby_idx1(1)).ID);
                        end
                        %Set c_range to mean +/- 1 standard deviation
                        c_range = [IntermodelStats(cfreq_idx).mean-IntermodelStats(cfreq_idx).std IntermodelStats(cfreq_idx).mean+IntermodelStats(cfreq_idx).std];
                        %Adjust to log/log10, if color data is log/log10
                        if strcmp(color_data_scale,'log'); c_range = log(c_range);
                        elseif strcmp(color_data_scale,'log10'); c_range = log10(c_range);
                        end
                    else
                        color_warn = true;
                        warning('LOG_PANEL:NoIntermStats',['No intermodel stats found for ',filevarFN,freq,...
                            ' ',expArray{colorby_idx1(1)},', reverting to default color range (from min to max of ',...
                            'color by dataset, determined individually by model'])
                        c_range = def_crange;
                    end
                else
                    c_range = def_crange;
                end
            end
            
            %If an element of the color range is unreal, switch back to
            %default
            if any(~isreal(c_range))
                color_warn = true;
                warning('LOG_PANEL:ImaginaryCRange',['Colorbar range contains nonreal values. ',...
                    'If no custom color range had been selected, this is likely due to non-normality ',...
                    'of color data, in which case using mean +/- 1 std ',...
                    'from IntermodelStats may give unphysical results, which can turn nonreal if ',...
                    'the color data is scaled logarithmically. Reverting to default color range (',...
                    'from min to max of color by dataset, determined individually by model.'])
                c_range = def_crange;
            end
            
            %Throw warning if incosistent colorscheme will be applied
            %across models
            if color_warn
                %For data_ratio and data, since this calculation is done
                %per model, will be inconsistent if just using the data
                %ranges themselves. Suggest manual set to avoid.
                warning('LOG_PANEL:MismatchedColors',...
                    ['Since no across-model manual color range has been identified, colormaps are only internally consistent within models. ',...
                    'The master colormap likely does not reflect the colormap used in every model individually, ',...
                    'since it is calibrated to the range of the last column. To avoid this, use the "color_range" flag and ',...
                    'manually set the color range across models.']);
            end
            
            %Change data scale, if desired, and adjust colorbar labels
            %accordingly; additionally, add latex mathtype markers ($)
            if strcmp(color_data_scale,'raw')
                def_clabel = ['$',def_clabel,'$']; %#ok<AGROW>
                %Do nothing, default behavior
            elseif strcmp(color_data_scale,'log')
                color_data = log(color_data);
                def_clabel = ['$\mathrm{log}\left(',def_clabel,'\right)$']; %#ok<AGROW>
            elseif strcmp(color_data_scale,'log10')
                color_data = log10(color_data);
                def_clabel = ['$\mathrm{log}_{10}\left(',def_clabel,'\right)$']; %#ok<AGROW>
            end
            
            %Throw warning if over 30% of points lie beyond color range
            if length(find(color_data<c_range(1) | color_data>c_range(2)))/numel(color_data) > 0.3
               warning('LOG_PANEL:ShortColorScale',[num2str(100*(length(find(color_data<c_range(1) | color_data>c_range(2)))/numel(color_data)),2),' percent ',...
                   'of colormap values lie beyond the color range [',num2str(c_range),']. If this is unintended, ',...
                   'adjust color range manually to a wider range.'])
            end
            
            %Set colors
            if isa(c_range,'numeric') && isequal(c_range,[0 1])
                cmap = flipud(autumn);
            elseif strcmp(color_by,'data_ratio') && ~isempty(strfind(color_data_scale,'log'))
                cmap = 'ColorRtGtB'; %Log ratio data goes Blue to Gray to Red
            elseif strcmp(color_by,'color')
                cmap = colorby_idx1;
            else
                cmap = 'def'; %Non-Ratio data uses a flipped jet scale
            end
            
            %Add units to clabel if data or log data is colored by
            if strcmp(color_by,'data')
                def_clabel = [def_clabel,' $',units,'$']; %#ok<AGROW>
            end
            
            %Set rgb color value for each point (x_data{vars},y_data{vars})
            [point_colors_tmp,cmap,c_range] = data_colormap(color_data,color_split,'range',c_range,'colormap',cmap);
            [point_colors{1:num_freqs}] = deal(point_colors_tmp); clear point_colors_tmp

            %Match indices of colormap with indices of data
            for vars = 1:num_freqs
                %Clip out points as with the data to graph above
                point_colors{vars}(idx_display_clip{vars},:) = [];
                %Set indices to be shuffled correctly with data
                point_colors{vars} = point_colors{vars}(shuffled_idxs{vars},:);
            end
            
            %Set overall figure colormap (same as subplots, used for master
            %colorbar)
            colormap(cmap)
            
            %% 4.3.2.6 Plot data by model / frequency band
            
            for vars = 1:num_freqs
                %% 4.3.1.6a Create Subplots
                if num_models > 1
                    if num_figures>1 %6 x num_freqs grid, for num_models > 6
                        subplot(num_freqs,6,model_vars_subscripts(vars,j))
                    else %num_models x num_freqs grid, for num_models <= 6
                        subplot(num_freqs,num_models,model_vars_subscripts(vars,j));
                    end
                elseif num_models == 1
                    subplot('Position',DisplaySettings(disp_idx).subp_position(vars,:))
                end
                
                hold on
            
                %% 4.3.1.6b Set Graph Dimensions
                ax = gca;
                axis(axis_limits)
                %Set axes ticks to the same; in integer steps
                ax.YTick = (axis_limits(3):axis_limits(4));
                ax.XTick = (axis_limits(1):axis_limits(2));
                ax.FontSize = DisplaySettings(disp_idx).axes_fontsize;
                %Make axes equal
                axis square
                
                %% 4.3.1.6c Plot Data
                %Scatterplot of data
                sc = scatter(x_data{vars}(:),y_data{vars}(:),point_size,point_colors{vars},'DisplayName',[num2str(clip_lat_min),'^o - ',num2str(clip_lat_max),'^o lat (absolute)']);
                sc.MarkerFaceColor = marker_face_color;
                %Area-weighted mean dot
                mudot = plot(x_data_mean{vars,(f-1)*6+j},y_data_mean{vars,(f-1)*6+j},'.k','MarkerSize',DisplaySettings(disp_idx).button_size,'DisplayName','mean (area-weighted)');
                %Print values of x and y mean (for funsies, mainly)
                fprintf(['x_mean = ',num2str(x_data_mean{vars,(f-1)*6+j}),'\ny_mean = ',num2str(y_data_mean{vars,(f-1)*6+j}),'\n']);
                
                %% 4.3.1.6d Principal Component Analysis
                idx_fit = cell(num_freqs,1);
                %Find points that are to be fit (not necessarily same as
                %displayed)
                idx_display = setdiff(1:length(lon)*length(lat),idx_display_clip{vars});
                idx_fit_tmp = setdiff(1:length(lon)*length(lat),idx_fit_clip{vars});
                %Get (unshuffled, unskipped) indices of points displayed
                %that shouldn't be fit through
                [~,idx_disp_fit_diff] = setdiff(idx_display,idx_fit_tmp);
                %Find out which of those indices are in the shuffled,
                %skipped (through idx_count) vector shuffled_idxs
                %(idx_fit_clip2 are the indices in shuffled_idxs that
                %correspond to idx_disp_fit_diff)
                [~,idx_fit_clip2] = ismember(idx_disp_fit_diff,shuffled_idxs{vars});
                %Set the final indices of x_data and y_data to be fit
                %through, based off of x_data/y_data being already
                %shuffled, clipped
                idx_fit{vars} = setdiff(1:length(x_data{vars}),idx_fit_clip2);
                %Store % of displayed points not included in fit for
                %subplot caption
                pct_fitclipped = (1-length(idx_fit{vars})/length(x_data{vars}))*100;
                
                %Perform PCA on points to be fit
                [coeff{vars,(f-1)*6+j},~,~,~,explained{vars,(f-1)*6+j},mu{vars,(f-1)*6+j}] = pca([x_data{vars}(idx_fit{vars});y_data{vars}(idx_fit{vars})]','Weights',area_weights{vars});

                %Get a linspace the width of plot
                X1 = linspace(axis_limits(1),axis_limits(2));
                %Calculate slope of 1st principal component
                slp{vars,(f-1)*6+j} = coeff{vars,(f-1)*6+j}(2,1)/coeff{vars,(f-1)*6+j}(1,1);
                
                %Plot 1st principal component
                if show_fit 
                    pca_line = plot(X1,slp{vars,(f-1)*6+j}*(X1-mu{vars,(f-1)*6+j}(1))+mu{vars,(f-1)*6+j}(2),'LineWidth',DisplaySettings(disp_idx).line_width,'DisplayName','1st principal component');
                    pca_line.Color = [128/255,141/255,169/255];
                end
                
                %% 4.3.1.6e Indicate off-scale points
                %Find off-scale idxs (beyond graph limits)
                outl{4} = find(y_data{vars} > axis_limits(4));
                outl{3} = find(y_data{vars} < axis_limits(3));
                outl{2} = find(x_data{vars} > axis_limits(2));
                outl{1} = find(x_data{vars} < axis_limits(1));
                %Set those values to the axis limits, so they are plotted
                %on the axis
                x_data{vars}(outl{2}) = axis_limits(2); x_data{vars}(outl{1}) = axis_limits(1);
                y_data{vars}(outl{4}) = axis_limits(4); y_data{vars}(outl{3}) = axis_limits(3);
                %Get all off-scale indices
                pts_offscale = unique([union(outl{4},outl{3}),union(outl{2},outl{1})]);
                %Get percentage of points offscale
                pct_offscale = 100*(numel(pts_offscale))/length(y_data{vars}(:));
                %Plot offscale points on axis 
                scatter(x_data{vars}(pts_offscale),y_data{vars}(pts_offscale),4,point_colors{vars}(pts_offscale,:),'filled');
                
                %% 4.3.1.6f Annotate Subplot
                %Graph reference y = x line
                line(X1,X1,'Color', 'k', 'LineStyle', '--','LineWidth',DisplaySettings(disp_idx).line_width)
                
                %Define subplot annotation
                Annotate_String = [];
                if show_fit; Annotate_String = [Annotate_String,'slope = ',num2str(round(slp{vars,(f-1)*6+j},2,'significant')),char(10)]; end %#ok<AGROW>
                Annotate_String = [Annotate_String,'corr = ',num2str(round(CorrCoeffs{vars,(f-1)*6+j}(2),2,'significant')),char(10)]; %#ok<AGROW>
                if pct_offscale > 0; Annotate_String = [Annotate_String,'offscale ',num2str(round(pct_offscale,2,'significant')),'%',char(10)];end %#ok<AGROW>
                if clip_minval && clip_minval_display; Annotate_String = [Annotate_String,'fit excludes ',num2str(round(pct_fitclipped,2,'significant')),'%',char(10)];end %#ok<AGROW>
                if clip_insig; Annotate_String = [Annotate_String,'outside confidence interval: ',num2str(round((num_ciclipped(vars)/area_points{vars})*100,2,'significant')),'%',char(10)];end %#ok<AGROW>
                
                %Insert annotation
                if annotate_subplot
                    annotation('textbox',ds2nfu(DisplaySettings(disp_idx).scaption_position),...
                        'String',Annotate_String,'FitBoxToText','on','FontSize',DisplaySettings(disp_idx).caption_fontsize,...
                        'LineStyle','none','VerticalAlignment','top','HorizontalAlignment','left');
                    clear Annotate_String
                end
                
                %Create legend
                if disp_idx ==4 && show_legend
                    if show_fit; lgd = legend([sc mudot pca_line]); 
                    else lgd = legend([sc mudot]);
                    end
                    lgd.Box = 'off';
                    lgd.Location = 'southeast';
                    lgd.FontSize = DisplaySettings(disp_idx).caption_fontsize;
                end
                
                %Create dotted cross at origin
                if origin_cross
                    line([origin_cross_dims(1),origin_cross_dims(2)],[0 0],'Color','k','LineStyle','--','LineWidth',DisplaySettings(disp_idx).line_width)
                    line([0 0],[origin_cross_dims(3),origin_cross_dims(4)],'Color','k','LineStyle','--','LineWidth',DisplaySettings(disp_idx).line_width)
                end
                
                %% 4.3.1.6g Set subplot axis labels, if desired
                %If plotting more than one model, title with model,
                %frequency descriptions
                if num_models>1
                    %Title columns with model names
                    if vars==1
                        try
                            model_disp = various_defaults.modelArray_disp(find(cellfun(@(x) strcmp(modelArray{j},x),various_defaults.modelArray_disp(:,1))),2); %#ok<FNDSB>
                        catch
                            model_disp = modelArray{j};
                        end
                        title(strcat(model_disp),'FontSize',DisplaySettings(disp_idx).title_fontsize);
                    end
                
                    %Title rows with frequency names
                    if j==1 && num_freqs > 1
                        freq_name = Data(struct_idxs(1)).Raw(freq_ids(1,vars)).Name; %extract name from frequency-dependent data struct (since matching, doesn't matter which one)
                        ylabel(freq_name,'FontWeight','bold','FontSize',DisplaySettings(disp_idx).title_fontsize); clear freq_name
                    end
                %If plotting one model / one freq, use axis labels
                elseif length(modelArray) == 1
                    xlabel(['\fontsize{',num2str(DisplaySettings(disp_idx).axes_fontsize),'}{0}',xlabel_string],'Interpreter','latex','FontSize',DisplaySettings(disp_idx).axes_fontsize);
                    ylabel(['\fontsize{',num2str(DisplaySettings(disp_idx).axes_fontsize),'}{0}',ylabel_string],'Interpreter','latex','FontSize',DisplaySettings(disp_idx).axes_fontsize);
                    if show_title
                        if num_freqs == 1
                            title([modelArray{j},' (',Data(struct_idxs(1)).Raw(freq_ids(1,vars)).Name,')'])
                        else
                            title(Data(struct_idxs(1)).Raw(freq_ids(1,vars)).Name);
                        end
                    end
                end
                
            end
            
        end
        
        %% 4.3.3 Create Master Colorbar
        if show_colorbar
            axes_colorbar = axes('Visible','off','Parent',figure1,...
                'Position',[0.05 0.05 0.95 0.95],...
                'Clim',c_range);
            c=colorbar(axes_colorbar,...
                'Position', DisplaySettings(disp_idx).colorbar_position);
            if isempty(cust_clabel)
                colorbarlabel = def_clabel;
            else
                colorbarlabel = cust_clabel;
            end
            c.FontSize = DisplaySettings(disp_idx).colorbar_fontsize;
            ylabel(c,colorbarlabel,'FontSize',DisplaySettings(disp_idx).colorbar_fontsize,'interpreter',cinterpreter);
        end
        
        %% 4.3.4 Annotate Figure
        %Create caption
        if show_caption
            if isempty(cust_caption)
                %Get long experiment names
                exp_name = cell(2,1);
                for i = 1:2
                    try
                        exp_name(i) = various_defaults.expArray_disp(find(cellfun(@(x) strcmp(expArray{i},x),various_defaults.expArray_disp(:,1))),2); %#ok<FNDSB>
                    catch %If not in various_defaults, just use the experiment 'raw'/save name
                        exp_name(i) = expArray(i);
                    end
                end
                %Set caption string
                Annotate_String = [freqdesc,vardesc,', ',exp_name{1},' vs. ',exp_name{2},' ',hemisphere,' ',season];
                if scale_byvar; Annotate_String = [Annotate_String,'; scaled by mean changes in ',freqdesc_scale,vardesc_scale]; end %#ok<AGROW>
                if region_clip; Annotate_String = [Annotate_String,char(10),'Region ',region_id]; end %#ok<AGROW>
                if lat_clip; Annotate_String = [Annotate_String,char(10),'Latitude ',num2str(clip_lat_min),'-',num2str(clip_lat_max),' (absolute)']; %#ok<AGROW>
                else Annotate_String = [Annotate_String,char(10), 'Latitude 0-90 (absolute)']; %#ok<AGROW>
                end
                if clip_insig; Annotate_String = [Annotate_String,char(10),'Data outside of 95\% confidence interval removed'];end %#ok<AGROW>
                if clip_minval; Annotate_String = [Annotate_String,'; removed below ',num2str(round(clip_minval_min,2,'significant')),' $',units,'$']; end %#ok<AGROW>
                if clip_maxval; Annotate_String = [Annotate_String,'; removed above ',num2str(round(clip_maxval_max,2,'significant')),' $',units,'$']; end %#ok<AGROW>
                if num_figures ~= 1; Annotate_String = [Annotate_String,char(10),'(Figure ',num2str(f),' of ',num2str(num_figures),')']; end %#ok<AGROW>
                if idx_count ~=1; Annotate_String = [Annotate_String,'; every ',num2str(idx_count),' point(s) (randomly chosen) plotted']; end %#ok<AGROW>
                if weight_pca; Annotate_String = [Annotate_String,'; PCA calculated with area-weighting of pixels']; end %#ok<AGROW>
            else
                Annotate_String = cust_caption;
            end
            %Print caption (setting fontsize within tex environment)
            annotation(figure1,'textbox',...
                        DisplaySettings(disp_idx).caption_position,...
                        'String',{['\fontsize{',num2str(DisplaySettings(disp_idx).annotate_fontsize),'}{0}',Annotate_String]},'Interpreter','latex',...
                        'FitBoxToText','off',...
                        'EdgeColor','none');
        end
        
        %Create cartoon axes to show axes labels
        if DisplaySettings(disp_idx).axes_cartoon
            axes_cartoon = axes('Parent',figure1,...
                'Position',[0.5 0.03 0.04 0.05],'FontSize',DisplaySettings(disp_idx).axes_cartoon_fontsize);
            xlabel(axes_cartoon,xlabel_string,'interpreter','latex')
            ylabel(axes_cartoon,ylabel_string,'interpreter','latex','Rotation',0,'HorizontalAlignment','right')
            axes_cartoon.XTick = [];
            axes_cartoon.YTick = [];
        end
        
        %% 4.3.5 Save Figure
        if save_graphs
            if f == 1
                %Add colormap identifier
                filename_base = strcat(filename_base,colorby_filename);
            end
            
            %Add figure number
            filename_fig = strcat(filename_dir,filename_base,'_',num2str(f),'of',num2str(num_figures));
            %Add additional filename suffixes (_TEST, etc.)
            filename_fig = [filename_fig,filename_add]; %#ok<AGROW>
            
            %Save graph as .eps
            print(gcf,'-depsc2','-loose',filename_fig);
            %Save graph as .png if desired
            if save_png
                print(gcf,'-dpng','-loose',filename_fig);
            end
            
            EndMsg=[filename_fig,' Complete'];
            disp(EndMsg);
        end
        
        %Close graph if desired
        if close_graphs
            close(figure1);
        end
        
    end
    
    %% 4.4 Save Stats
    if save_stats
        modelArray = [modelArray_sub{:}];
        if ~exist([filename_dir,'Stats/'])
            mkdir([filename_dir,'Stats/'])
            warning([filename_dir,'Stats/ created.'])
        end
        filename_stats = strcat(filename_dir,'Stats/',filename_base,'_stats',filename_add);
        if ~isempty(struct_idxs); freq_names{:} = {Data(struct_idxs(1)).Raw(freq_ids(1,:)).Name};
        else freq_names = 'no frequency-dependent data involved';end 
        save(filename_stats,'-v7.3','modelArray','freq_names','CorrCoeffs','slp','mu','y_data_mean','x_data_mean','coeff','explained','xlabel_string','ylabel_string')
        EndMsg=[filename_stats,' Complete'];
        disp(EndMsg);
    end
    
end

%% 5 Cleanup
%Restore supressed warnings
warning('on','DATA_COLORMAP:TooFewColors')


