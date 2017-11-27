function varargout = clip_var_min(min_toclip,data_set,varargin)
% CLIP_VAR_MIN  get indices below a minimum value.
%   idx_minval = CLIP_VAR_MIN(min_toclip,data_set) returns the indices of 
%   the pixels of [data_set] with a value below [min_toclip] - behavior is
%   identical to idx_minval = find(data_set < min_toclip).
%
%   [idx,idx_minval,num_totclipped,num_clipped] = 
%   CLIP_VAR_MIN(min_toclip,data_set,idx) adds the indices below
%   [min_toclip] to existing index group for more complex clipping. 
%   All indices not in [idx_minval] are added to [idx], with [num_clipped] 
%   (the number of unique indices added to [idx]) added to
%   [num_totclipped]. [idx_minval] is outputted as a cell array of
%   dimensions length(idx) x 1.
%
%   [min_toclip] should be entered in units of the data (units are
%   described in [code_dir]/Varnames.csv or through the function
%   var_chars).
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also CLIP_VAR_MAX, CLIP_LAT, CLIP_REGION, CLIP_CI, FIND
%
%   For questions/comments, contact Kevin Schwarzwald
%   (kschwarzwald@uchicago.edu)
%   Last edit: 08/25/2016

if (~isempty(varargin))
    idx = varargin{1};
    num_vars = length(idx);
    idx_minval = cell(num_vars,1);
    num_totclipped = zeros(num_vars,1);
    num_clipped = zeros(num_vars,1);
    
    for vars = 1:num_vars
        idx_minval{vars} = find(data_set < min_toclip);
        if (~isempty(varargin))
            pre_length = length(idx{vars});
            idx{vars} = unique(cat(1,idx{vars},idx_minval{vars}));
            num_totclipped(vars) = length(idx{vars});
            num_clipped(vars) = length(idx{vars}) - pre_length;
        end
    end
else
    %idx_minval = find(reshape(data_set, [numel(data_set) 1]) < min_toclip);
    idx_minval = find(data_set < min_toclip);
end

if nargout == 1
    varargout{1} = idx_minval;
else
    if isempty(varargin)
        error('More than one output requires an input of [idx]. Not enough inputs.')
    end
    varargout{1} = idx;
    varargout{2} = idx_minval;
    varargout{3} = num_totclipped;
    varargout{4} = num_clipped;
end
