function varargout = clip_var_max(max_toclip,data_set,varargin)
% CLIP_VAR_MAX  get indices above a maximum value.
%   idx_maxval = CLIP_VAR_MAX(max_toclip,data_set) returns the indices of 
%   the pixels of [data_set] with a value above [max_toclip] - behavior is
%   identical to idx_maxval = find(data_set > max_toclip).
%
%   [idx,idx_maxval,num_totclipped,num_clipped] = 
%   CLIP_VAR_MAX(max_toclip,data_set,idx) adds the indices below
%   [max_toclip] to existing index group for more complex clipping. 
%   All indices not in [idx_maxval] are added to [idx], with [num_clipped] 
%   (the number of unique indices added to [idx]) added to
%   [num_totclipped]. [idx_maxval] is outputted as a cell array of
%   dimensions length(idx) x 1.
%
%   [max_toclip] should be entered in units of the data (units are
%   described in [code_dir]/Varnames.csv or through the function
%   var_chars).
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also CLIP_VAR_MIN, CLIP_LAT, CLIP_REGION, CLIP_CI, FIND
%
%   For questions/comments, contact Kevin Schwarzwald
%   (kschwarzwald@uchicago.edu)
%   Last edit: 08/25/2016

if (~isempty(varargin))
    idx = varargin{1};
    num_vars = length(idx);
    idx_maxval = cell(num_vars,1);
    num_totclipped = zeros(num_vars,1);
    num_clipped = zeros(num_vars,1);
    
    for vars = 1:num_vars
        idx_maxval{vars} = find(data_set > max_toclip);
        if (~isempty(varargin))
            pre_length = length(idx{vars});
            idx{vars} = unique(cat(1,idx{vars},idx_maxval{vars}));
            num_totclipped(vars) = length(idx{vars});
            num_clipped(vars) = length(idx{vars}) - pre_length;
        end
    end
else
    %idx_maxval = find(reshape(data_set, [numel(data_set) 1]) < max_toclip);
    idx_maxval = find(data_set > max_toclip);
end

if nargout == 1
    varargout{1} = idx_maxval;
else
    if isempty(varargin)
        error('More than one output requires an input of [idx]. Not enough inputs.')
    end
    varargout{1} = idx;
    varargout{2} = idx_maxval;
    varargout{3} = num_totclipped;
    varargout{4} = num_clipped;
end
