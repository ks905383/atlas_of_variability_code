function varargout = clip_ci(var_idx,model,expArray,varargin)
% CLIP_CI  isolate indices without meaningful magnitude change.
%   
%   idx_ci = CLIP_CI(var_idx,model,expArray) returns the indices of the
%   pixels for [model] and variable identifier [var_idx] determined to
%   HAVE meaningful changes in variability between two experiment runs
%   given by [expArray] (2-sigma, as calculated in Variability_CI).
%
%   [idx,idx_ci,num_totclipped,num_clipped] =
%   CLIP_CI(var_idx,model,expArray,idx) adds the insignificant pixels to
%   existing index group for more complex clipping. All indices deemed NOT
%   meaningful (this is the opposite behavior as [idx_ci] above) are added
%   to [idx], with [num_clipped] (the number of pixels in [idx_ci] uniquely
%   added to [idx], i.e. how many of [idx_ci] were not already in [idx])
%   added to [num_totclipped].
%
%   Function requirements: _StdDevCI.mat file output from Variability_CI;
%   functions structfind (from MathWorks community), var_chars, name_chars.
%
%   CLIP_CI(...,'[flag]',[params],...) modify program run as below:
%       'use_freqs',[cell]      - get indices corresponding to the
%                                 frequency bins described in [cell]. These 
%                                 can either be in ID form
%                                 {'HF','HFavgDev',etc.} or frequency bin
%                                 size {3,15,90,etc.}.
%       'filename_add',[str]    - search for a file with extra non-standard
%                                 identifiers (i.e. '_TEST'); pairs to
%                                 [filename_add] in Variability_CI. The
%                                 [str] is added after '_StdDevCI' to the
%                                 filename to load.
%
%   This function is paired with Variability_CI in that it primarily
%   interprets output from that function saved in '_StdDevCI' files. As
%   such, IDs and frequency bin sizes need to be equal to those generated
%   in Variability_CI. 
%
%   NOTE: this function is part of the Atlas of Variability code package
%
%   All directories listed as [____] are set in various_defaults.m
%
%   See also: VARIABILITY_CI, CLIP_LAT, CLIP_REGION
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last edit: 10/05/2016

%Get rid of some name_chars warnings
warning('off','name_chars:MultConvs'); %Unnecessary if you're clear on what you're doing

%Manage inputs
various_defaults = matfile('various_defaults.mat');
filename_add = [];
spec_freqs = false;
if nargout > 1
    idx = varargin{1};
    varargin{1} = 0; %To avoid switch tripping over cell array
end
if (length(varargin)>1) %all of these take multiple inputs, so if one, it's [idx]
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'filename_add'}
                filename_add = varargin{in_idx+1};
            case {'use_freqs'}
                spec_freqs = true;
                use_freqs = varargin{in_idx+1}; 
                varargin{in_idx+1} = 0; %To avoid switch tripping over cell array
        end
    end
end

%Get file characteristics
[~,filevar,freq,~,~,~,~,~,~] = var_chars(var_idx);
[~,run,~,~,~] = name_chars(model,expArray,filevar,freq);

%Load _StdDevCI file
filename_c = strcat(various_defaults.proc_data_dir,model,'/',filevar,freq,model,'_',expArray{1},'_',expArray{2},'_',run{1},'StdDevCI',filename_add,'.mat');
load(filename_c);

if ~spec_freqs %if using all frequency bins in given _StdDevCI file
    if nargout > 1 && length(StdDevsRatio) ~= length(idx)
        errMsg = ['Significance calculated for different number of frequencies as currently ',...
                    'used in analysis (taken from length of [idx] cell input); cannot provide accurate clipped indices (to fix, ',...
                	'make sure the _StdDevCI file contains data for the same number of ',...
            		'frequencies as desired in the analysis or choose specific frequencies with [use_freqs] flag).'];
        error('clip_ci:InexistFreq',errMsg)
    end
    %Set struct indices to all indices in struct
    id_idxs = (1:length(StdDevsRatio))';
else %if using specific frequencies within _StdDevCI file
   %Support both frequency bin size and frequency bin id inputs
   if isa(use_freqs{1},'char')
       struct_searchfield = 'ID';
   elseif isa(use_freqs{1},'numeric')
       struct_searchfield = 'Size';
   else
       error('clip_ci:InvalidIds','use_freqs must be a cell array of ID strings {HF, LF, etc.} or a cell array of frequency bin sizes {3, 15, 90, etc.}');
   end
   %Pre-allocate struct indices
   id_idxs = zeros(length(use_freqs),1);
   for vars = 1:length(use_freqs)
       try
           id_idxs(vars) = structfind(StdDevsRatio,struct_searchfield,use_freqs{vars});
       catch
           if isa(use_freqs{vars},'numeric'); use_freqs{vars} = num2str(use_freqs{vars}); end
           error('clip_ci:InexistFreq',['File does not have data for frequency bin given by [',use_freqs{vars},']']);
       end
   end
end
       
idx_ci_tmp = cell(length(id_idxs),1);
for vars = 1:length(id_idxs)
    %Get indices of non-significant pixels (excuse the eval call) 
    idx_ci_tmp{vars} = sub2ind(size(StdDevsRatio(id_idxs(vars)).Data), StdDevsRatio(id_idxs(vars)).lon_ciidx,StdDevsRatio(id_idxs(vars)).lat_ciidx);
end

if nargout == 1
    idx_ci = cell(length(id_idxs),1);
    for vars = 1:length(id_idxs)
        full_idxs = (1:numel(StdDevsRatio(vars).Data))';
        %Get indices of meaningful pixels, if only one output
        idx_ci{vars} = full_idxs(~ismember(full_idxs,idx_ci_tmp{vars}));
    end
    varargout{1} = idx_ci;
    if mod(length(varargin),2) == 1;
        warning('clip_ci:UnnecClip','idx was inputted, but only idx_ci will be outputted (showing meaningful pixels) due to only one output requested');
    end
end

if nargout > 1 && (~isempty(varargin))
    %Collate with existing subsets
    num_totclipped = zeros(length(idx),1);
    num_clipped = zeros(length(idx),1);
    num_ciclipped = zeros(length(idx),1);
    idx_ci = idx_ci_tmp;
    for vars = 1:length(id_idxs)
        %Get number of clipped before the adding of lat clipped indices
        pre_length = length(idx{vars});
        %Add new indices to be removed to idx
        idx{vars} = unique(cat(1,idx{vars},idx_ci{vars}));
        num_totclipped(vars) = length(idx{vars});
        num_clipped(vars) = length(idx{vars}) - pre_length;
        num_ciclipped(vars) = length(idx_ci);
    end
    varargout{2} = idx_ci;
    varargout{1} = idx;
    varargout{3} = num_totclipped;
    varargout{4} = num_clipped;
end