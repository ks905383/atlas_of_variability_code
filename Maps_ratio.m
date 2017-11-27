function Maps_ratio(VarIndices,modelArray,varargin)
% MAPS_RATIO    Map changes in variability.
%
%   MAPS_RATIO(VarIndices,modelArray) creates a set of maps showing the
%   ratio of variability between two experiments (by default RCP8.5 and
%   piControl) globally of variables given by [VarIndices] for models given
%   by [modelArray], separated by frequency bands. By default, all
%   frequency bands but data for 'all freqs' (the basic standard deviation)
%   are plotted. Maps are saved as .eps files in
%   [figure_dir]/FinalFigs/[VarDomain]/.
%
%   Resulting figure is by default configured for 11" x 8" (792 x 576 pts).
%
%   Save and data locations and experiment names are taken from
%   various_defaults.mat
%
%   Required data: a *SpectralRatios.mat file with data in the
%   [proc_data_dir] from [various_defaults.m] containing the ratios of the
%   desired variables; if using polygon mapping, then a .nc file in the
%   [raw_data_dir] with the same lat/lon grid as the data, containing
%   variables [lat_bnds] and [lon_bnds] (if these are not available,
%   lat/lon bounds are constructed from prevailing lat/lon grid, details
%   below).
%
%   MAPS_RATIO(...,'[flag]',[params],...) modify program run as below:
%       'experiments',[cell]- set experiments between which ratio is
%                             plotted (requires _sqrtPower for both
%                             experiments) (def; {'rcp85','piControl'})
%       'show_freqs',[cell] - set frequency bands to graph, either 'all' or
%                             a cell array of either desired frequency band
%                             sizes ([cell] = {5,30,365...}) or IDs ([cell]
%                             = {'HF','XF',...}). (def: set to graph all
%                             but full (all frequencies) bands). If only
%                             one frequency is used, the filename is
%                             extended by the '_[ID]'.
%       'mark_ci'           - mark pixels not meaningfully different from 1
%                             with a gray cross (additionally requires a
%                             _StdDevCI file for each variable/model)
%       'ci_marker_size',[int] - set size of ci cross (default 0.2)
%       'seasons',[cell]    - set seasonal data to be plotted (i.e.
%                             {'DJF'})
%       'scale_by',[char]   - by default 'log10', set the colormap scaling.
%                             Ticks will still be in raw units, but the
%                             scale will change. Options are 'log10',
%                             'log', and 'raw'. 
%       'close_graphs'      - close graphs in UI after saving
%       'visible',['off'/'on']    - set whether to display figure in UI
%                                   (def 'on')
%       'pointsy',[int]     - change resultant image size (height)
%       'pointsx',[int]     - change resultant image size (width)
%       'save_png',[log] - set whether to save .png copy using gs
%                             (def: false)
%       'save_fig',[log]    - set whether to save .fig copy (def: false)
%       'replace',[log]     - set whether to replace existing files (def:
%                             true)
%       'clim',[arr],([arr])- set colorbar limits (def. [-0.5 0.5] if
%                             scaled by log ([0.6 1.6] in original units)
%                             or log10 (default, [0.36 3.6] in original
%                             units), or [0.5 1.5] if not scaled). If a
%                             second array is inputted, it is interpreted
%                             as the location of the colorbar ticks. By
%                             default, tickmarks are set at [0.01 0.02 0.05
%                             0.1 0.2 0.5 1 2 5 10 20 50 100] (only [0.5 1
%                             2] are visible in the default log10, clim =
%                             [-0.5 0.5] setting), or [0 0.5 1 1.5 2] for
%                             scale_by = 'raw'). 
%       'color_split',[log],[int] - set whether to split colormap into
%                                   discrete colors (default: true, 32)
%       'testing'           - adds '_TEST' to filename
%       'save_log'          - exports log of attempted saves with
%                             success/failure status under
%       'save_dir',[char]   - set directory to save outputted figures in,
%                             by default set to the output directory
%                             specified by various_defaults. 
%                             ~/CalcLogs/Maps_pi_flat_log_[starttime].txt
%       'lim_map',[log/int] - input map lat/lon limits. If set to 'true',
%                             then the max and min lat/lon of the data are
%                             used. If you intend to manually set limits,
%                             instead input two 2-element vectors, one for
%                             lat limits, one for lon limits (i.e., in the
%                             default case, 'lim_map',[-90 90],[-180 180]).
%                             (setting to 'false' shows a world map).
%                             Filename is changed to include the limits.
%       'plot_method',[char]    - set method of plotting on map. By
%                                 default, this code uses the matlab
%                                 function 'pcolorm', which is fast but
%                                 removes a row and column from each side
%                                 (pixels are plotted with positions as
%                                 vertices, see for example <a href =
%                                 "http://www.mathworks.com/matlabcentral/fileexchange/50706-offsets-and-missing-data-via-pcolor-and-surf?focused=3873981&tab=example">here</a>).
%                                 The other option is to input 'polygon',
%                                 which plots every pixel as an individual
%                                 polygon. This method fixes the pcolor
%                                 distortion issue, but is much slower,
%                                 requires extra data, and currently
%                                 features a possibly solved bug (NOT FIXED
%                                 FOR LAT/LON 2D ARRAYS!) where for some
%                                 LongRunMIP models, the resultant figure
%                                 is entirely NaN-colored (the bug is
%                                 related to the direction of the lat
%                                 vector). The extra data necessary is a
%                                 .nc file in the raw data folder with the
%                                 same number of lat/lons as the data that
%                                 also contains the variables lat_bnds and
%                                 lon_bnds. If this file does not exist,
%                                 lat and lon bounds are created by simply
%                                 assuming every pixel spans halfway
%                                 between the relevant coordinate and
%                                 adjacent coordinates in all directions
%                                 (in which case a warning will be thrown).
%                                 Since this method assumes that if the
%                                 length of the lat_bnds variable is equal
%                                 to the length of the lat variable (same
%                                 with lon), issues could arise if for some
%                                 reason there exist two separate grids
%                                 with the same number of elements for the
%                                 same model. This method also implicitly
%                                 assumes compliance with the CMIP5
%                                 defaults (rectangular pixels in lat/lon
%                                 space with edges saved as
%                                 'lat/lon_bnds'), which is not always the
%                                 case for non-CMIP5 variables.
%       'break_on_error',[log]  - if true, then errors will cause program
%                                 run to end. By default = false, and
%                                 errors are caught, with warnings
%                                 outputted, but run continues.
%       'transpose_data',[log]  - set whether to transpose array before
%                                 plotting. By default, true. It's a bit
%                                 unclear why this is necessary, but it's a
%                                 hack that makes sure the resultant map is
%                                 oriented correctly. In very rare cases,
%                                 the resultant map(s) might be flipped, in
%                                 which case, rerun with
%                                 'transpose_data',false to fix.
%                                 
%
%
%   Saving convention:
%   [filevar]_[freq]_[model]_[exp1]_[exp2]_StdDevRatio_Maps...
%                                                       (_CImarked)...
%                                           (_[lat lims]_[lon lims)...
%                                                                     .eps
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   KNOWN BUGS/ETC.: If there are weird smears, white lines, etc. in your
%   map, it's probably related to the process of wrapping to 360. In
%   general, pcolorm has issues at the wraparound break (the last/first
%   lon/data column), often leaving a white line. A fix exists generally
%   for when this gap is around the Prime Meridian, but issues still exist
%   when that is not the case. Currently includes a hacky solution. Please
%   see section 3.2.3.
%
%   See also
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified 06/26/2017

%% 1. Set Defaults
various_defaults = matfile('various_defaults.mat');
%Set clock (for logging purposes)
format shortG
startTimestamp = clock;

%Set data defaults
expArray = {'rcp85','piControl'};
seasArray = {'all'};
use_freqs = {'all_but_full'};
all_freqs = false;
season_folder = various_defaults.season_dir;

%Set figure defaults
plot_method = 'pcolor';
lim_map = false; map_lat_lim = []; map_lon_lim = [];
mark_ci = false;
ci_marker_size = 0.2;
pointsy = 576;
pointsx = 792;
close_graphs = false;
visible = 'on';
color_split = true;
clim_tmp = [-0.5 0.5];
Ticks = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100];
scale_by = 'log10';
clabel = [];
colorbar_bands = 32;
sort = true;
show_caption = true;

%Set mapping defaults
map_origin = 155;
projection = 'bsam';
transpose_data = true;

%Set size defaults
annotation_fontsize = 10;
colorbar_fontsize = 16;
colorbar_ticks_fontsize = 12;
title_fontsize = 11;

%Set saving defaults
save_log = false;
save_eps = true;
save_png = false;
save_fig = false;
replace_files = true;
filename_add = [];
alt_save_dir = [];

%Set run defaults
break_on_error = false;

%Set behavior of optional function flags
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'seasons'}
                seasArray = varargin{in_idx+1};
                varargin{in_idx+1} = 0;
            case {'experiments'}
                expArray = varargin{in_idx+1};
                varargin{in_idx+1} = 0; %To stop switch from tripping over a cell
            case {'show_freqs'}
                use_freqs = varargin{in_idx+1};
                if isa(use_freqs,'char') && strcmp(use_freqs,'all')
                    all_freqs = true;
                end
                varargin{in_idx+1} = 0; %To stop switch from tripping over a cell
                sort = false; %To make sure it prints in order inputted
            case {'scale_by'}
                scale_by = varargin{in_idx+1};
            case {'replace'}
                replace_files = varargin{in_idx+1};
            case {'close_graphs'}
                close_graphs = true;
            case {'mark_ci'}
                mark_ci = true;
            case {'clabel'}
                clabel = varargin{in_idx+1};
            case {'clim'}
                clim_tmp = varargin{in_idx+1}; varargin{in_idx+1} = 0;
                if isa(varargin{in_idx+2},'numeric')
                    Ticks = varargin{in_idx+2}; varargin{in_idx+2} = 0;
                end
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
            case {'save_png'}
                save_png = varargin{in_idx+1};
            case {'save_fig'}
                save_fig = varargin{in_idx+1};
            case {'save_eps'}
                save_eps = varargin{in_idx+1};
            case {'testing'}
                filename_add = [filename_add,'_TEST']; %#ok<AGROW>
            case {'color_split'}
                color_split = varargin{in_idx+1};
                colorbar_bands = varargin{in_idx+2};
                if isa(varargin{in_idx+2},'numeric')==0 || varargin{in_idx+2} < 0
                    error('Invalid colorbar split, must be a positive integer')
                end
            case {'visible'}
                visible = varargin{in_idx+1};
            case {'filename_add'}
                filename_add = [filename_add,varargin{in_idx+1}]; %#ok<AGROW>
            case {'plot_method'} %Options = pcolor, polygon
                plot_method = varargin{in_idx+1};
            case {'lim_map'}
                if isa(varargin{in_idx+1},'logical')
                    lim_map = varargin{in_idx+1};
                elseif isa(varargin{in_idx+1},'numeric')
                    lim_map = true;
                    map_lat_lim = varargin{in_idx+1}; varargin{in_idx+1} = 0;
                    map_lon_lim = varargin{in_idx+2}; varargin{in_idx+2} = 0;
                end
            case {'map_origin'}
                map_origin = varargin{in_idx+1};
            case {'projection'}
                projection = varargin{in_idx+1};
            case {'save_dir'}
                alt_save_dir = varargin{in_idx+1};
            case {'ci_marker_size'}
                ci_marker_size = varargin{in_idx+1};
            case {'colorbar_fontsize'}
                colorbar_fontsize = varargin{in_idx+1};
            case {'colorbar_ticks_fontsize'}
                colorbar_ticks_fontsize = varargin{in_idx+1};
            case {'title_fontsize'}
                title_fontsize = varargin{in_idx+1};
            case {'show_caption'}
                show_caption = varargin{in_idx+1};
            case {'break_on_error'}
                break_on_error = varargin{in_idx+1};
            case {'transpose_data'}
                transpose_data = varargin{in_idx+1};
        end
    end
end

%Make use_freqs, modelArray a cell if it's not for code below to work with 
%general type inputs
if ~iscell(use_freqs)
    use_freqs = {use_freqs};
end
if ~iscell(modelArray)
    modelArray = {modelArray};
end

%% 2. Setup Graph Characteristics

%Set experiment markers
exp_marker = cell(2,1);
for i = 1:2
    try
        exp_marker(i) = various_defaults.expArray_disp(find(cellfun(@(x) strcmp(expArray{i},x),various_defaults.expArray_disp(:,1))),2); %#ok<FNDSB>
    catch
        exp_marker(i) = expArray(i);
    end
end

%Define colormap
standard_colormaps = matfile('standard_colormaps');
%Set colormap
if color_split
    colormap_ratios = standard_colormaps.ColorRtWtB(1:(256/colorbar_bands):256,:); 
    %Make sure colormap is white on both sides of 0
    colormap_ratios = [colormap_ratios(2:colorbar_bands/2+1,:); colormap_ratios(colorbar_bands/2+1,:); colormap_ratios(colorbar_bands/2+2:end,:)];
else
    colormap_ratios = standard_colormaps.ColorRtWtB; 
end

%Set colorbar and tickmarks
if strcmp(scale_by,'log10')
    lTicks = log10(Ticks); c_lim = clim_tmp;
elseif strcmp(scale_by,'log')
    lTicks = log(Ticks); c_lim = clim_tmp;
elseif strcmp(scale_by,'raw')
    if isequal(Ticks,[0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100]); lTicks = [0 0.5 1 1.5 2]; else; lTicks = Ticks; end
    if isequal(clim_tmp,[-0.5 0.5]); c_lim = [0.5 1.5]; else; c_lim = clim_tmp; end 
end
%Round tickmarks to 3 sig figs
Ticks = round(Ticks,3,'significant');

%Pre-allocate log array
if save_log
    complete_log = cell(1,length(VarIndices)*length(modelArray),length(seasArray));
end

%Get number of models
num_models = length(modelArray);

%% 3. Plot PI Diagnostics Maps
for i = 1:length(VarIndices)
    %Get variable characteristics
    [~,filevarFN,freq,vardesc,~,~,~,freqdesc,~] = var_chars(VarIndices(i));
    for j=1:num_models
        for seas = 1:length(seasArray)
            initial_vars = who;
            try
                %% 3.1 Model-Specific Graph Setup
                %Get run characteristics from filenames
                if strcmp(seasArray{seas},'all')
                    season_file_add = []; folder = '/'; season_desc = [];
                else; season_file_add = ['_',seasArray{seas}]; folder = ['/',season_folder,'/'];
                    season_desc = seasArray{seas};
                end
                
                filename_t = [various_defaults.proc_data_dir,modelArray{j},folder,filevarFN,freq,modelArray{j},'_',expArray{1},'_',expArray{2},'_','*_SpectralRatios',season_file_add,'.mat'];
                file_t = dir(filename_t);
                if isempty(file_t); error('MapsStdDev:NoFile',['File ',filename_t,' not found, no data available']); end
                filename = file_t.name;
                delim = strcat('\_','|','\.'); splitStr = regexp(filename,delim,'split');
                run{1} = strcat(splitStr{7},'_');
                [~,run] = name_chars(modelArray{j},expArray,filevarFN,freq);
                
                %Set figure output directory
                if ~isempty(alt_save_dir)
                    save_dir = alt_save_dir;
                else
                    save_dir = [various_defaults.figure_dir,FileDomain(filevarFN),'/'];
                end
                filenameS = [save_dir,filevarFN,freq,modelArray{j},'_',expArray{1},'_',expArray{2},'_',run{1},'StdDevRatio_Maps',season_file_add];
                if mark_ci
                    filenameS = strcat(filenameS,'_CImarked');
                end
                if ~isempty(map_lat_lim)
                    filenameS = strcat(filenameS,'_',num2str(map_lat_lim),'_',num2str(map_lon_lim));
                end
                
                %% 3.2 Create maps if desired
                if replace_files || (~replace_files && ~exist([filenameS,'.eps'],'file'))
                    %Load Ratios
                    load([various_defaults.proc_data_dir,modelArray{j},folder,filename]);
                    
                    %Fix polygon bug with different latitude vector
                    %conventions that causes entire figure to be grey
                    if strcmp(plot_method,'polygon')
                       if lat(1) > 0
                           lat = flipud(lat);
                           for idx = 1:length(StdDevsRatio)
                              StdDevsRatio(idx).Data = fliplr(StdDevsRatio(idx).Data);  %#ok<AGROW>
                           end
                       end
                    end
                    
                    %% 3.2.1 Get lat/lon bounds, if using polygons to graph
                    if strcmp(plot_method,'polygon')
                        nc_files = dir([various_defaults.raw_data_dir,modelArray{j},'/*.nc']);
                        derive_bnds = true;
                        if ~isempty(nc_files)
                            n = 1; file_nlat = []; file_nlon = [];
                            %Cycle through candidate .nc files to find
                            %matching lat/lon grid (by length of lats/lons)
                            while ~(isequal(file_nlat,size(lat,1)) && isequal(file_nlon,size(lon,1)))
                                %Break if reached past number of possible
                                %.nc files to scan
                                if n > length(nc_files)
                                    %Since derive_bnds is already true, the
                                    %warning below will take care of it
                                    break
                                end
                                
                                %Look at a candidate .nc file, determine if
                                %the file contains lat_bnds and lon_bnds,
                                %and see if nlat/nlon match
                                ncfile_info = ncinfo([various_defaults.raw_data_dir,modelArray{j},'/',nc_files(n).name]);
                                if ~all(cellfun(@(x) ~strcmp('lat_bnds',x),{ncfile_info.Variables.Name}))
                                    %Get number of lat/lon values in nc file
                                    file_latbnds_size = ncfile_info.Variables(cellfun(@(x) strcmp('lat_bnds',x),{ncfile_info.Variables.Name}) ).Size;
                                    file_nlat = file_latbnds_size(file_latbnds_size~=2);
                                    file_lonbnds_size = ncfile_info.Variables(cellfun(@(x) strcmp('lon_bnds',x),{ncfile_info.Variables.Name}) ).Size;
                                    file_nlon = file_lonbnds_size(file_lonbnds_size~=2);
                                end
                                
                                %If nlat/nlon match with loaded variables,
                                %load lat and lon bounds
                                if isequal(file_nlat,size(lat,1)) && isequal(file_nlon,size(lon,1))
                                    lat_bnds = ncread([various_defaults.raw_data_dir,modelArray{j},'/',nc_files(n).name],'lat_bnds');
                                    lon_bnds = ncread([various_defaults.raw_data_dir,modelArray{j},'/',nc_files(n).name],'lon_bnds');
                                    %Note that lat/lon bounds don't have to
                                    %be derived
                                    derive_bnds = false;
                                    %Make sure individual pixels are rows
                                    %instead of columns
                                    if find(file_latbnds_size~=2) == 2
                                        lat_bnds = lat_bnds';
                                    end
                                    if find(file_lonbnds_size~=2) == 2
                                        lon_bnds = lon_bnds';
                                    end
                                end
                                %Increase while-loop counter
                                n = n+1;
                            end
                        end
                        
                        %If need to interpolate lat/lon bounds, do so
                        if derive_bnds
                            warning('MapsStdDev:NoBnds',['Latitude/longitude bounds for ',modelArray{j},...
                                ' that match lat/lon grid for ',filevarFN,freq,...
                                ' were not found in ',various_defaults.raw_data_dir,...
                                '. Lat/lon grid edges will be interpolated from loaded lat/lon grid.'])
                            %Pre-allocate lat/lon bound array
                            lat_bnds = zeros(length(lat),2);
                            lon_bnds = zeros(length(lon),2);
                            %Have support for different lat/lon counters
                            if lat(1) < 0; lat_bnds(1,1) = -90; else; lat_bnds(1,1) = 90; end
                            %Interpolate lat boundaries by taking the
                            %average of subsequent pixels
                            for lat_idx = 1:length(lat)-1
                                lat_bnds(lat_idx,2) = (lat(lat_idx+1)-lat(lat_idx))/2+lat(lat_idx);
                                if lat_idx < length(lat) && lat_idx>1
                                    lat_bnds(lat_idx,1) = lat_bnds(lat_idx-1,2);
                                end
                            end
                            %Set last element
                            lat_bnds(length(lat),1) = lat_bnds(length(lat)-1,2);
                            if lat(length(lat)) < 0; lat_bnds(length(lat),2) = -90; else; lat_bnds(length(lat),2) = 90; end
                            
                            %Have support for different lat/lon counters
                            if lon(1) < 0; lon_bnds(1,1) = -180; elseif lon(1) ==0; lon_bnds(1,1) = 0; elseif lon(1)>0; lon_bnds(1,1) = 180;end
                            %Interpolate lat boundaries by taking the
                            %average of subsequent pixels
                            for lon_idx = 1:length(lon)-1
                                lon_bnds(lon_idx,2) = (lon(lon_idx+1)-lon(lon_idx))/2+lon(lon_idx);
                                if lon_idx < length(lon) && lon_idx>1
                                    lon_bnds(lon_idx,1) = lon_bnds(lon_idx-1,2);
                                end
                            end
                            %Set last element
                            lon_bnds(length(lon),1) = lon_bnds(length(lon)-1,2);
                            if lon(length(lon)) < 0; lon_bnds(length(lon),2) = -180;
                            elseif lon(length(lon)) == 0; lon_bnds(length(lon),2) = 0;
                            elseif lon(length(lon)) < 180; lon_bnds(length(lon),2) = 180;
                            elseif lon(length(lon)) > 180; lon_bnds(length(lon),2) = 360;
                            end
                        end
                        
                        
                        %Generate struct for each pixel
                        clear struct_poly
                        struct_poly(length(lon)*length(lat)).Lat = 0;  %#ok<AGROW>
                        for poly_idx = 1:length(lon)*length(lat)
                            [lon_idx,lat_idx] = ind2sub([length(lon) length(lat)],poly_idx);
                            struct_poly(poly_idx).Lat = [lat_bnds(lat_idx,1); lat_bnds(lat_idx,2); lat_bnds(lat_idx,2); lat_bnds(lat_idx,1); lat_bnds(lat_idx,1);NaN];
                            struct_poly(poly_idx).Lon = [lon_bnds(lon_idx,1); lon_bnds(lon_idx,1); lon_bnds(lon_idx,2); lon_bnds(lon_idx,2); lon_bnds(lon_idx,1);NaN];
                            struct_poly(poly_idx).Geometry = 'Polygon';
                        end
                    end
                    
                    %% 3.2.2. Set up frequency bands/statistics to graph
                    %Get frequency band characteristics
                    num_vars = length(StdDevsRatio); %#ok<*NODEF>
                    graph_titles = cell(num_vars,1); freq_band_sizes = zeros(num_vars,2); freq_band_ids = cell(num_vars,1);
                    for band = 1:num_vars
                        graph_titles{band} = StdDevsRatio(band).Name;
                        if ~isempty(StdDevsRatio(band).Size)
                            freq_band_sizes(band,:) = StdDevsRatio(band).Size;
                        else
                            freq_band_sizes(band,:) = [0 0];
                        end
                        freq_band_ids{band} = StdDevsRatio(band).ID;
                    end
                    
                    %Decide which frequency bands to graph
                    if ~all_freqs
                        if isa(use_freqs{1},'char') && strcmp(use_freqs{1},'all_but_full')
                            if contains([freq_band_ids{:}],'Full')
                                [~,bands,~] = setxor(freq_band_ids,{'Full'});
                            else
                                bands = 1:length(StdDevsRatio);
                            end
                        elseif isa(use_freqs{1},'char')
                            for s_idx = 1:length(use_freqs)
                                if ~isempty(structfind(StdDevsRatio,'ID',use_freqs{s_idx}))
                                    bands(s_idx) = structfind(StdDevsRatio,'ID',use_freqs{s_idx});
                                else; error('MapsStdDev:NoFreq',['The frequency band identified by ',use_freqs{s_idx},' does not exist in the SpectralRatios file.',...
                                        ' Please choose one or more of the following IDs: ',strjoin({StdDevsRatio(:).ID})])
                                end
                            end
                        elseif isa(use_freqs{1},'numeric')
                            for s_idx = 1:length(use_freqs)
                                bands(s_idx) = structfind(StdDevsRatio,'Size',intersect(freq_band_sizes,cell2mat(use_freqs(s_idx,:))));
                            end
                        end
                        if ~strcmp(use_freqs{1},'all_but_full')
                            filenameS = [filenameS,'_',use_freqs{:}]; %#ok<AGROW>
                        end
                    else
                        bands = 1:length(StdDevsRatio);
                    end
                    
                    %Sort bands to use in ascending length by frequency band
                    %size
                    if sort
                        try
                            [~,sort_idxs] = sort(cell2mat({StdDevsRatio(bands).Size}'),1);
                            bands = bands(sort_idxs); clear sort_idxs
                        catch
                            warning('MapStdDev:FailedSort','Sorting subplots by frequency band size unsuccessful. Will print in order of struct.')
                        end
                    end
                        
                    
                    %Load ci pixel indices, if desired
                    if mark_ci
                        [ci_idxs,~,~,~] = clip_ci(VarIndices(i),modelArray{j},expArray,cell(length(bands),1),'use_freqs',use_freqs);
                        [lonxx,latyy] = ndgrid(lon,lat);
                    end
                    
                    %% 3.2.3 Plot Data
                    %Create figure
                    figure1 = figure('visible',visible,'Colormap',colormap_ratios);
                    %Get World Map
                    coast = load('coast');
                    
                    %Get subplot positions based on number of frequencies
                    %graphed
                    position_set = subplot_pos(length(bands));
                    
                    %Wrap data around 360 degrees to avoid white line from
                    %pcolorm's origin behavior (works together with
                    %wrapping data around below in plot_data = [plot_data'
                    %plot_data(1,:)']');
                    if size(lon,2) == 1 
                        if lon(end)>345
                            lon = [lon' 360'];
                        else 
                        %(with hack to make sure crazy things don't happen
                        %if lon doesn't start counting from 0). Ideally
                        %this should also prevent the white line, but this
                        %doesn't seem to fully be the case yet...
                            if lon(1) < 0; first_last(1) = 360+lon(1); else; first_last(1) = lon(1); end
                            if lon(end) < 0; first_last(2) = 360+lon(end); else; first_last(2) = lon(end); end
                            lon = [lon' mean(first_last)];
                        end
                        lonWrapped = wrapTo360(lon);
                    else
                        lonWrapped = lon;
                    end
                    
                    %Plot data by frequency band / subplot
                    for band = 1:length(bands)
                        %Take Logs
                        if strcmp(scale_by,'log10')
                            plot_data = log10(StdDevsRatio(bands(band)).Data);
                        elseif strcmp(scale_by,'log')
                            plot_data = log(StdDevsRatio(bands(band)).Data);
                        elseif strcmp(scale_by,'raw')
                            plot_data = StdDevsRatio(bands(band)).Data;
                        end
                        
                        if strcmp(plot_method,'pcolor')
                            %Wrap data around 360 degrees
                            plot_data = [plot_data' plot_data(1,:)']';
                            
                            if transpose_data
                                plot_data = plot_data.';
                            end
                        elseif strcmp(plot_method,'polygon')
                            %Get color attributes
                            plot_colors = data_colormap(plot_data,colorbar_bands,'colormap',colormap_ratios,'range',c_lim);
                            color_atts = makesymbolspec('Polygon',{'INDEX',[1 numel(struct_poly)],'FaceColor',plot_colors,'LineStyle','none'});
                        else
                            error('MapsRatio:IncorrectPlotMethod','Incorrect plot method. Must be "polygon" or "pcolor"')
                        end
                        
                        
                        %Set up subplot
                        subplot('position',position_set(band,:),'Visible','off','Parent',figure1,...
                            'FontSize',8,'XTick',[],'YTick',[],'DataAspectRatio',[1 1 1],...
                            'clim',c_lim);
                        %Get geographic map axes with set projection, origin
                        if ~lim_map
                            axesm(projection,'Origin',map_origin)
                        else
                            if isempty(map_lat_lim)
                                map_lat_lim = [double(min(lat(:))) double(max(lat(:)))];
                                map_lon_lim = [double(min(lon(:))) double(max(lon(:)))];
                                [~,tmp_inda] = min(lon(:));
                                [~,tmp_indb] = max(lon(:));
                                if tmp_indb<tmp_inda
                                    map_lon_lim = fliplr(map_lon_lim);
                                    map_lon_lim(map_lon_lim<0) = 360-(180+map_lon_lim(map_lon_lim<0));
                                end
                            end
                            axesm(projection,'MapLatLimit',...
                                map_lat_lim,'MapLonLimit',map_lon_lim);
                        end
                        %Plot data onto map
                        if strcmp(plot_method,'pcolor')
                            pcolorm(double(lat),double(lonWrapped),plot_data);
                        elseif strcmp(plot_method,'polygon')
                            geoshow(struct_poly,'SymbolSpec',color_atts)
                        end
                        %Plot coasts on top of map
                        geoshow(coast.lat,coast.long,'Color','black')
                        framem off; gridm off; axis off; mlabel off; plabel off;
                        %Title
                        title(graph_titles{bands(band)},'FontSize',title_fontsize);
                        %Mark unmeaningful pixels, if desired
                        if mark_ci && ~isempty([latyy(ci_idxs{band}) lonxx(ci_idxs{band})])
                            plotm([latyy(ci_idxs{band}) lonxx(ci_idxs{band})],'+','MarkerSize',ci_marker_size,'MarkerEdgeColor',[0.2 0.2 0.2])
                        end
                    end
                    
                    %% 3.2.4 Create Master Colorbar
                    axes1 = axes('Visible','off','Parent',gcf,...
                        'Position',[0.05 0.05 0.95 0.95],...
                        'CLim',c_lim);
                    if length(bands)==1 %Set to align with map 
                        cposition = [0.91 0.15-(0.2/14.8)*0.75 0.014 0.75*(14.3/14.8)];
                    else
                        cposition = [0.91 0.15 0.014 0.75];
                    end
                    c=colorbar('peer',axes1,...
                        'Position',cposition);
                    set(c,'Ytick',lTicks,'YTicklabel',Ticks,'fontsize',colorbar_ticks_fontsize);
                    if isempty(clabel) 
                        if length(use_freqs) == 1 && strcmp(use_freqs{1},'mean')
                            colorbarlabel = ['$\mu_{',exp_marker{1},'}/\mu_{',exp_marker{2},'}$'];
                        elseif length(use_freqs) > 1 && any(cellfun(@(x) strcmp(x,'mean'), use_freqs))
                            colorbarlabel = ['$',exp_marker{1},'/',exp_marker{2},'$'];
                        else
                            colorbarlabel=['$\sigma_{',exp_marker{1},'}/\sigma_{',exp_marker{2},'}$'];
                        end
                    else %If custom colorbar label
                        colorbarlabel = clabel;
                    end
                    if VarIndices(i) == 31
                        colorbarlabel = ['$\frac{1}{\Delta \mu}\frac{\sigma_{',exp_marker{1},'}}{\sigma_{',exp_marker{2},'}}$'];
                    end
                    ylabel(c,colorbarlabel,'interpreter','latex','FontSize',colorbar_fontsize);
                    
                    %% 3.2.5 Annotate Figure
                    if show_caption
                        %Set caption descriptions
                        AnnotateString = [freqdesc,vardesc,' ',season_desc,newline,modelArray{j},newline];
                        if mark_ci
                            AnnotateString = [AnnotateString,'Gray crosses show pixels not meaningfully different from 1']; %#ok<AGROW>
                        end
                        
                        %Set annotation
                        annotation(figure1,'textbox',...
                            [0.05 0.01 0.9 0.0710172744721689],...
                            'String',{AnnotateString},...
                            'FitBoxToText','off',...
                            'FontSize',annotation_fontsize,...
                            'EdgeColor','none');
                    end
                    
                    %% 3.2.6 Print Figure
                    %Fix bug in print that causes white polygons at edge of
                    %plot to turn black upon printing
                    if strcmp(plot_method,'polygon')
                       figure1.InvertHardcopy = 'off';
                       figure1.Color = [1 1 1];
                    end
                    
                    figure1.PaperUnits = 'points';
                    set(gcf,'PaperPositionMode','manual');
                    figure1.PaperPosition = [0 0 pointsx pointsy];
                    if ~isempty(filename_add)
                        filenameS = [filenameS,filename_add]; %#ok<AGROW>
                    end
                    
                    if ~exist(save_dir,'file') || save_eps || save_fig || save_png
                        mkdir(save_dir)
                        warning([save_dir,' was not found, and was created.'])
                    end
                    
                    if save_eps
                        print(gcf,'-depsc2','-loose',filenameS);
                    end
                    disp([filenameS,' saved!'])
                    if save_fig
                        saveas(figure1,filenameS)
                    end
                    if save_png
                        print(gcf,'-dpng','-loose',filenameS);
                    end
                    EndMsg=['Std Dev Ratio Maps for ',modelArray{j},' ',seasArray{seas},' ',filevarFN,freq,' Complete'];
                    disp(EndMsg);
                    
                    %Store success message in log
                    if save_log
                        save_log_msg = [filevarFN,freq,modelArray{j},seasArray{seas},' Complete'];
                        complete_log{(i-1)*length(modelArray)+j,seas} = save_log_msg;
                    end
                    
                else
                    disp(['File for ',modelArray{j},' ',filevarFN,freq,seasArray{seas},' already exists, was not replaced']);
                    %Store message that no figure needed in log
                    if save_log
                        save_log_msg = [filevarFN,freq,modelArray{j},' already exists'];
                        complete_log{(i-1)*length(modelArray)+j,seas} = save_log_msg;
                    end
                end
            catch ME
                disp(ME)
                disp(ME.message); disp(ME.stack); disp([ME.stack.line])
                %Store warning in log
                if save_log
                    save_log_wrn = ['****',filevarFN,freq,modelArray{j},seasArray{seas},' Not Complete****'];
                    complete_log{(i-1)*length(modelArray)+j,seas} = save_log_wrn;
                end
                if break_on_error
                    error('MAPS_RATIO:break',['Graph for ',seasArray{seas},' ',modelArray{j},' ',filevarFN,freq,' not created!'])
                else
                    warningMsg = ['Graph for ',seasArray{seas},' ',modelArray{j},' ',filevarFN,freq,' not created!'];
                    warning(warningMsg);
                end
            end
            if close_graphs
                close(gcf)
            end
            clearvars('-except',initial_vars{:})
        end
    end
    EndMsg=['Std Dev Ratio Maps for ',filevarFN,freq,' for All Models Complete'];
    disp(EndMsg);
end
%% 4. Export log as table
if save_log
    complete_log = complete_log';
    complete_log = complete_log(:);
    export_log = cell2table(complete_log','VariableNames',{'VariabilityMaps_attempted_executions'});
    writetable(export_log,['/home/kschwarzwald/CalcLogs/MAPS_RATIO_log_',num2str(startTimestamp(1)),'-',num2str(startTimestamp(2)),'-',num2str(startTimestamp(3)),'-',num2str(startTimestamp(4)),'-',num2str(startTimestamp(5)),'-',num2str(startTimestamp(6)),'.txt']);
end
end
