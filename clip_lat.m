function varargout = clip_lat(lat_min,lat_max,lat,lon,varargin)
% CLIP_LAT  isolate indices in a certain latitude band.
%   idx_bnd = CLIP_LAT(lat_min,lat_max,lat,lon) returns the indices of the
%   pixels in the (absolute) latitude band given by [lat_min] and [lat_max]
%   in the resolution given by [lat] and [lon]
%
%   [idx_bnd,subset_params] = CLIP_LAT(lat_min,lat_max,lat,lon)
%   also outputs the exact lat and long bounds in a struct subset_params
%   (since resolutions differ between models, the closest latitude bands to
%   lat_min and lat_max are chosen; the exact value of those latitudes are
%   thus saved).
%
%   [idx,subset_params,idx_bnd,num_totclipped,num_clipped] =
%   CLIP_LAT(lat_min,lat_max,lat,lon,idx) adds all indices NOT part of
%   [lat_min,lat_max] band to existing index group for more complex
%   clipping. All indices not in [idx_bnd] (this [idx_bnd] is the opposite,
%   indexing-wise, than the one in above code syntaxes; represents pixels
%   NOT in lat band) are added to [idx], with [num_clipped] (the length of
%   [idx_bnd]) added to [num_totclipped].
%
%   By default, absolute latitude bands are outputted (so, lat_min and
%   lat_max values of 55,90 will give indices for [-90,-55] and [55,90]).
%   To only output indices for one hemisphere, 
%   [___] = CLIP_LAT(___,'hemisphere',['NH'/'SH']) with any of the above
%   syntaxes will do so.
%
%   [lat] and [lon] are nlat/nlon x 1 vectors, as taken from CMIP5 raw
%   netcdf files. [idx] is a num_vars x 1 cell array of indices to be added
%   to.
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   See also CLIP_REGION
%
%   For questions/comments, contact Kevin Schwarzwald
%   (kschwarzwald@uchicago.edu)
%   Last edit: 09/12/2016
one_sided = false;
if (~isempty(varargin))
    if length(varargin)==1 || length(varargin)==3
        idx = varargin{1};
    end
    if length(varargin)>2
        varargin = varargin(2:end);
    end
    if strcmp(varargin{1},'hemisphere')
        one_sided = true;
        h = varargin{2};
    %else
     %   error('clip_lat:InvalidInput',[varargin{1},' is not a supported option flag']);
    end
    varargout = cell(5,1);
else
    if nargout == 2
        varargout = cell(2,1);
    elseif nargout <= 1
        varargout = cell(1,1);
    end
end

nlat = length(lat);
nlon = length(lon);

%Remove latitude bands from poles and equator (if desired)
subset_params = struct;

%% Identify lat band closest to desired equatorial clip
if lat_min ~= 0;
    
    %This code assumes that latitudes are counted S-N. If N-S, flip. This
    %should not affect output if absolute latitude is asked for. 
    if lat(1)>0 && ~one_sided
        lat = flipud(lat);
    end
    
    tmp_eq_max = abs(lat-lat_min);
    tmp_eq_min = abs(lat+lat_min);
    
    [~,n_bnd] = min(tmp_eq_max); clear tmp_eq_max
    [~,s_bnd] = min(tmp_eq_min);
    %Quick hack to maintain symmetric clipping if multiple
    %minima
    if length(tmp_eq_min)>s_bnd && tmp_eq_min(s_bnd) == tmp_eq_min(s_bnd+1);
        s_bnd = s_bnd+1; clear tmp_eq_min
    end
    
    %Populate longitude coordinates
    lon_bnds_tmp = repmat((1:nlon),(n_bnd-s_bnd+1),1);
    lon_bnds_eq = lon_bnds_tmp(:); clear lon_bnds_tmp
    
    %Populate latitude coordinates, counting from equator
    lat_bnds_eq = repmat((s_bnd:n_bnd),1,nlon);
    subset_params.eq_bnds = {length(s_bnd:n_bnd)};
    subset_params.eq_lats = {lat(s_bnd),lat(n_bnd)};
    %Calculate indices of clipped bands from equator
    idx_bnds_eq = sub2ind([nlon nlat],lon_bnds_eq,lat_bnds_eq');
    
else %If no clip at equator desired
    idx_bnds_eq = [];
    subset_params.eq_bnds = {[]};
    subset_params.eq_lats = {[]};
end

%% Identify lat band closest to desired polar clip
if lat_max ~= 90 && lat_max ~= max(lat)+(lat(2)-lat(1))/2;
    
    tmp_pl = abs(lat-lat_max);
    [~,clip_lat_pl_tmp] = min(tmp_pl); clear tmp_pl
    if clip_lat_pl_tmp >= (nlat/2 + 1)
        clip_lat_pl = nlat - clip_lat_pl_tmp; clear clip_lat_pl_tmp
    else %Make sure that clip_lat_pl is less than half the length of the lat vector, to allow for latitude vectors to be either increasing or decreasing
        clip_lat_pl = clip_lat_pl_tmp; clear clip_lat_pl_tmp
    end
    
    %Populate latitude coordinates, couting from poles (if/else accounts
    %for the fact that sometimes latitudes are saved descending...)  
    lat_bnds_pl = repmat([1:clip_lat_pl (nlat-clip_lat_pl+1):nlat],1,nlon);
    
    %Populate longitude coordinates
    lon_bnds_tmp = repmat((1:nlon),2*clip_lat_pl,1);
    lon_bnds_pl = reshape(lon_bnds_tmp,numel(lon_bnds_tmp),1);
    
    %Calculate indices of clipped bands from poles
    idx_bnds_pl = sub2ind([nlon nlat],lon_bnds_pl,lat_bnds_pl');
    subset_params.pl_bnds = {clip_lat_pl};
    subset_params.pl_lats = {lat(clip_lat_pl),lat(end-clip_lat_pl+1)};
    
else %If no clip at poles desired
    idx_bnds_pl = [];
    subset_params.pl_bnds = {[]};
    subset_params.pl_lats = {[]};
end

%% Create unified index vector
idx_bnds = cat(1,idx_bnds_pl,idx_bnds_eq);
%Get indices of points actually in the lat band
idx_bnd = find(~ismember((1:nlon*nlat)',idx_bnds));

if one_sided
    [~, latyy] = ndgrid(lon', lat');
    if strcmp(h,'NH')
        idx_bnds = find(latyy<0);
        idx_bnd = find(latyy>0);
    else
        idx_bnds = find(latyy>0);
        idx_bnd = find(latyy<0);
    end
end

%Set up data for outputs
varargout{1} = idx_bnd;
if nargout > 1
    varargout{2} = subset_params;
end

if exist('idx','var')
    %Collate with existing subsets
    num_totclipped = zeros(length(idx),1);
    num_clipped = zeros(length(idx),1);
    num_latclipped = zeros(length(idx),1);
    for vars = 1:length(idx)
        %Get number of clipped before the adding of lat clipped indices
        pre_length = length(idx{vars});
        %Add new indices to be removed to idx
        idx{vars} = unique(cat(1,idx{vars},idx_bnds));
        num_totclipped(vars) = length(idx{vars});
        num_clipped(vars) = length(idx{vars}) - pre_length;
        num_latclipped(vars) = length(idx_bnds);
    end
    varargout{3} = idx_bnds;
    varargout{1} = idx;
    varargout{4} = num_totclipped;
    varargout{5} = num_clipped;
end

end
