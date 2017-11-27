function varargout = clip_region(region_id,varargin)
% CLIP_REGION   get indices of a geographical region.
%   idx_region = CLIP_REGION(region_id,model,lat,lon,...) returns the
%   indices of the region given by [region_id] (input options described
%   below) in resolution of [model] given by [lat] and [lon]. [lat] and
%   [lon] can either be nlat x 1 and nlon x 1 vectors (basic rectangular
%   grid, with a lat,lon value for each row, column of data) or
%   equal-dimension nlon x nlat arrays (more complex geographic grid, with
%   a lat,lon pair at each point of data).
%
%   [idx_region,subset_params] = CLIP_REGION(region_id,model,lat,lon,..)
%   also outputs the exact lat and lon bounds of the original region
%   rectangle (before land/ocean masking) in a struct subset_params (since
%   resolutions differ between models, the closest lat/lon values to the
%   region boundaries are chosen; the exact value of those points are thus
%   saved). If [lat] and [lon] are vectors, each field of subset_params is
%   a 1x2 vector containing (1) the exact lat/lon used and (2) the index of
%   that value. If [lat] and [lon] are arrays, each field is a [num
%   selected lats/lons]x2 array containing the same.
%
%   [idx_region,subset_params,region_id] =
%   CLIP_REGION(region_id,model,lat,lon,...) also outputs the string region
%   identified [region_id] (a three-letter abbreviation of the region name)
%
%   [idx,subset_params,idx_region,num_totclipped,num_clipped,region_id] =
%   CLIP_REGION(region_id,model,lat,lon,idx,..)
%   adds all indices not part of region to existing index group for more
%   complex clipping. All indices not in [idx_region] are added to [idx],
%   with [num_clipped] (the length of [idx_region]) added to
%   num_totclipped. [region_id] is an optional add as above.
%
%   region_id = CLIP_REGION(region_idx) returns the string [region_id] given
%   by the input [region_idx], which can be a numerical region index
%   (giving a region's position in the struct contained in 'regions.mat',
%   the base file for this function).
%
%   region_name = CLIP_REGION(region_idx,'long') returns the long name of
%   the region in [region_name] (if present in struct), same input as
%   above.
%
%   [region_id] can either be a numerical index (of the struct Regions, found
%   in /Code/regions.mat), the three-letter id of the desired region, taken
%   from Castruccio et al. 2014 (Statistical Emulation of Climate Model
%   Predictions...), or a cell array of the above (i.e.
%   {'NEU','SEU','MED'}). [model] is the name of the desired model, as used
%   in the filing system /project/moyer/. [lat] and [lon] are nlat/nlon x 1
%   vectors, as taken from CMIP5 raw netcdf files. [idx] is a num_vars x 1
%   cell array of indices to be added to.
%   Some notable regions (apologies for the Eur/NA-centeredness):
%       Numerical Index     ID      Description
%                     1     G       Global (all points, no clip)
%                     2     GL      Global Land
%                     3     GO      Global Ocean
%                    11     WNA     West N America
%                    12     CNA     East C America
%                    13     ENA     East N America
%                    17     NEU     Northern Europe
%                    18     SEU     Southern Europe
%                    36     EPW     East Pacific West
%
%   Option (...,'size_z',size_z) gives instead the indices of a 3D-array
%   with dimensions [nlon nlat size_z] (with subsetting still done only in
%   lat and lon dimensions). This is useful if subsetting a global time
%   series (with a time dimension of size size_z) or a dataset with a
%   vertical dimension, etc.
%
%   Functions required (on path): subset_find; structfind (from community)
%
%   Data requirements: the file 'regions.mat'; if at least one of the
%   desired regions is of type 'OCEAN' or 'LAND', a file 'sftlf*.mat' in
%   the folder [proc_data_dir]/[model]/ with a nlon x nlat array 'sftlf'
%   giving the percentage of each pixel that is land for that model
%   (CLIP_REGION now supports multiple grids per model, as long as the
%   corresponding land fraction file is of the form "sftlf*.mat" in the
%   folder CLIP_REGION will cycle through files until the grid matches the
%   [lon x lat] inputted - this inputted [lon x lat] grid can cover a
%   smaller area as the 'sftlf*.mat' files, as long as the grid is the same
%   ('sftlf*.mat' files also include the lat, lon grid they are associated
%   with, to be able to compare and 'find' the smaller grid in the
%   larger)).
%
%   To add regions: edit the struct 'Regions' contained in
%   [code_dir]/regions.mat. Input max and min lat and lon values in the
%   corresponding fields, and, in the field 'Type', whether the region is
%   'LAND', 'OCEAN', or 'ALL' (if 'all', no contours are followed, the
%   resulting region is just a lat-lon rectangle). 'ID' is an easily
%   searchable string (such as 'SAH' for Sahara).
%
%   WARNING: Slight deviations from 'correct' geographical boundaries may
%   occur depending on model resolutions (i.e. a pixel in the South
%   Atlantic (SAT) showing up in a region subset of the East South Pacific
%   (SPE) if the resolution is low to enough to cause longitudinal grid
%   overlap between the coast off Tacna, Peru and the coast off Rio
%   Gallegos, Argentina). 
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also CLIP_LAT, CLIP_CI, CLIP_VAR_MIN,
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified: 01/31/2017

%Load region definitions and model-specific land/ocean mask
load('/project/moyer/Kevin/Code/regions.mat')

if isempty(varargin) || (length(varargin)==1 && strcmp(varargin{1},'long')) %If just getting region id/name string
    %Get index of desired region (can be either direct index input - numeric -
    %or region ID - string)
    if isa(region_id,'numeric')
        region_idx = region_id;
    elseif isa(region_id,'char')
        region_idx = structfind(Regions,'ID',region_id);
    end
    if isempty(varargin)
        if ~isa(region_id,'cell')
            varargout{1} = Regions(region_idx).ID;
        else
            varargout{1} = [];
        end
    else
        if ~isempty(Regions(region_idx).Name)
            varargout{1} = Regions(region_idx).Name;
        else
            warning('Long name not present for inputted [region_idx]. ID returned instead.')
            varargout{1} = Regions(region_idx).ID;
        end
    end
    
else %If getting full indices
    %% Setup
    if ~isa(region_id,'cell')
        region_id = {region_id};
    end
    
    model = varargin{1};
    lat = varargin{2};
    lon = varargin{3};
    varargin_opts = varargin(4:end);
    
    %Get optional function flags
    if (~isempty(varargin_opts))
        if length(varargin_opts)==1
            idx = varargin_opts{1};
        elseif length(varargin_opts)==2 && strcmp(varargin_opts{1},'size_z')
            size_z = varargin_opts{2};
        elseif length(varargin_opts)==3
            if strcmp(varargin_opts{1},'size_z')
                size_z = varargin_opts{2};
                idx = varargin_opts{3};
            elseif strcmp(varargin_opts{2},'size_z')
                size_z = varargin_opts{3};
                idx = varargin_opts{1};
            end
        else
            error('Wrong number of inputs. 4-6 needed.')
        end
        varargout = cell(6,1);
    else
        if nargout > 1 && nargout < 4
            varargout = cell(nargout,1);
        elseif nargout <= 1
            varargout = cell(1,1);
        end
    end
    
    %Get model resolution
    if size(lat,2) == 1 %If vector coordinates
        nlon = length(lon);
        nlat = length(lat);
    else %If array coordinates
        nlon = size(lat,1);
        nlat = size(lat,2);
    end
    
    %Turn W / E longitude into single, east-counting coordinates (to
    %support both 0 to 360 and -180 to 180 longitude standards)
    if any(lon(:) < 0)
        lon(lon < 0) = 360 + lon(lon < 0);
    end
    
    %% Load model-specific land/ocean mask
    %Figure out if a land/ocean mask is even needed (if all desired regions
    %are of type 'ALL' (just in a lat/lon bounding box), then the mask,
    %sftlf, is not needed)
    reg_types = cell(length(region_id),1);
    for reg = 1:length(region_id)
        if isa(region_id{reg},'numeric')
            region_idx = region_id{reg};
        elseif isa(region_id{reg},'char')
            region_idx = structfind(Regions,'ID',region_id{reg});
        end
        reg_types{reg} = Regions(region_idx).Type;
    end
    if all(cellfun(@(x) strcmp(x,'ALL'),reg_types))
        need_sftlf = false;
    else
        need_sftlf = true;
    end
    
    if need_sftlf
        sftlf_file = matfile(['/project/moyer/Kevin/',model,'/sftlf.mat']);
        sftlf = sftlf_file.sftlf;
        
        %Trim landmask to correct grid, if not already at the same grid
        if ~isequal(size(sftlf),[nlon nlat])
            %Get all sftlf files, if more than one
            search_files = dir(['/project/moyer/Kevin/',model,'/sftlf*.mat']);
            n=1;
            %Cycle through sftlf files till one is found with the correct grid
            while ~isequal(size(sftlf),[nlon nlat])
                if n>length(search_files)
                    error('CLIP_REGION:UnmatchedGrid',['The inputted lat/lon ',...
                        'arrays do not match the grid of any land fraction ',...
                        'variable found in the folder /project/moyer/Kevin/',...
                        model,'/. As a result, the region clip would be incorrect.',...
                        ' Please create a sftlf*.mat file with a sftlf array ',...
                        'with the same grid spacing as the one you are using right now ',...
                        '(as determined by the lat/lon vectors/arrays) before proceding.'])
                end
                sftlf_file = matfile(['/project/moyer/Kevin/',model,'/',search_files(n).name]);
                
                [sftlf_lat,sftlf_lon] = subset_find(sftlf_file.lat,sftlf_file.lon,lat,lon);
                sftlf = sftlf_file.sftlf(sftlf_lon,sftlf_lat);
            end
        end
    end
    
    %% Get indices of each region separately
    idx_region = cell(length(region_id),1);
    for reg = 1:length(region_id)
        %% Set region
        if isa(region_id{reg},'numeric')
            region_idx = region_id{reg};
        elseif isa(region_id{reg},'char')
            region_idx = structfind(Regions,'ID',region_id{reg});
        end
        
        %Start subset_params struct;
        subset_params.desc(1,:) = {'degree','index'};
        
        %% Find region indices
        %Get max/min lat/lon values wanted for given region
        coor_range = [Regions(region_idx).lat_min,Regions(region_idx).lat_max;...
            Regions(region_idx).lon_min,Regions(region_idx).lon_max];
        
        %Get indices of lat/lon values in model closest to those defined in Regions
        if size(lat,2) == 1 %If lat/lon are vectors
            [~,coor_lat_min] = min(abs(lat-coor_range(1,1)));
            [~,coor_lat_max] = min(abs(lat-coor_range(1,2)));
            [~,coor_lon_min] = min(abs(lon-coor_range(2,1)));
            [~,coor_lon_max] = min(abs(lon-coor_range(2,2)));
            
            %Store exact lat/lon used for subsetting (before land/ocean mask)
            if nargout > 1 && nargout < 4
                subset_params.lon_min(reg,:) = [lon(coor_lon_min),coor_lon_min];
                subset_params.lon_max(reg,:) = [lon(coor_lon_max),coor_lon_max];
                subset_params.lat_min(reg,:) = [lat(coor_lat_min),coor_lat_min];
                subset_params.lat_max(reg,:) = [lat(coor_lat_max),coor_lat_max];
            end
            
            %If lons cross prime meridian, make sure indices wrap around the right way
            if coor_lon_max < coor_lon_min
                lons = [coor_lon_min:size(lon) 1:coor_lon_max]';
            else
                lons = [coor_lon_min:coor_lon_max]'; %#ok<NBRAK>
            end
            
            %Get linear indices of the points in the rectangle of lat/lons found above
            idxs = sub2ind([nlon nlat],repmat(lons,1,length((coor_lat_min:coor_lat_max))),repmat((coor_lat_min:coor_lat_max),length(lons),1));
            idxs = idxs(:);
        else %If lat/lon are arrays
            lat_idxs = find(lat > coor_range(1,1) & lat < coor_range(1,2));
            lon_idxs = find(lon > coor_range(2,1) & lon < coor_range(2,2));
            idxs = intersect(lat_idxs,lon_idxs);
            
            %Store exact lat/lon used for subsetting (before land/ocean mask)
            if nargout > 1 && nargout < 4
                subset_params.lons(reg,:) = [lon(lon_idxs) lon_idxs];
                subset_params.lats(reg,:) = [lat(lat_idxs) lat_idxs];
            end
        end
        
        %Find land and ocean indices (defined as >50 land, <=50 ocean)
        if need_sftlf;
            land_idxs = find(sftlf>50);
            ocean_idxs = find(sftlf<=50);
        end
        
        %Remove land or ocean tiles for ocean or land type region, respectively
        if strcmp(Regions(region_idx).Type,'LAND')
            idx_region{reg} = intersect(idxs,land_idxs);
        elseif strcmp(Regions(region_idx).Type,'OCEAN')
            idx_region{reg} = intersect(idxs,ocean_idxs);
        elseif strcmp(Regions(region_idx).Type,'ALL')
            idx_region{reg} = idxs;
        else
            error(['"',Regions(region_idx).Type,'" not supported as a region type. Must be "LAND", "OCEAN", or "ALL"'])
        end
        
        %If 3d array, repeat idx subset for length of z
        if exist('size_z','var')
            idx_region_tmp = idx_region{reg};
            for i = 1:size_z-1
                idx_region{reg} = cat(1,idx_region{reg},(idx_region_tmp+i*(nlon*nlat)));
            end
        end
    end
    
    %% Combine all region indices into one master index vector
    idx_region = unique(cat(1,idx_region{:}));
    
    %% Output data
    %Set up data for outputs
    varargout{1} = idx_region;
    if nargout > 1 && nargout < 4
        varargout{2} = subset_params;
        if nargout == 3
            varargout{3} = Regions(region_idx).ID;
        end
    end
    
    if exist('idx','var')
        %Collate with existing subsets
        num_totclipped = zeros(length(idx),1);
        num_clipped = zeros(length(idx),1);
        for vars = 1:length(idx)
            idx_clip = find(~ismember((1:nlon*nlat)',idx_region));
            %Add new indices to be removed to idx
            idx{vars} = unique(cat(1,idx{vars},idx_clip));
            num_totclipped(vars) = length(idx{vars});
            num_clipped(vars) = length(idx_region);
        end
        varargout{3} = idx_region;
        varargout{1} = idx;
        varargout{4} = num_totclipped;
        varargout{5} = num_clipped;
        varargout{6} = Regions(region_idx).ID;
    end
    
end
end