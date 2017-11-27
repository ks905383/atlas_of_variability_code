function varargout = data_colormap(dataset,color_split,varargin)
% DATA_COLORMAP  get RGB values for a split (non-continuous) colormap for
%                every element in a dataset.
%
%   [rgb_adj] = DATA_COLORMAP(dataset,color_split) returns a
%   numel(dataset)x3 RGB vector [rgb_adj] containing a colormap based on
%   'jet' flipped up-down, split by the number of colors given by
%   color_split' for [dataset]. [rgb_adj] can be used as a colormap input.
%   Literally, the function rounds all values in [dataset] *DOWN* to the
%   closest (ceil(max(dataset))-floor(min(dataset))/color_split values, and
%   maps them to the colormap (so [-5,0)==-5, [5,10) == 5,etc.).
%
%   [rgb_adj,cmap,c_range] = DATA_COLORMAP(dataset,color_split) also
%   returns the unique rgb values used as the colormap in [cmap] and the
%   range of values represented by the colormap in a 2x1 array [c_range].
%
%   This function gives slightly more control over the splitting process
%   than simply [colormap](color_split) --> it is primarily intended to
%   easily map data to a colorbar for scatterplotting purposes.
%
%   Treating NaNs: By default, NaNs are set to a light grey [0.9 0.9 0.9].
%   Run methods can be changed below in options. 
%
%   Function requirements (on path): rgb (from MathWorks community, if
%   inputting CSS3 color names), piecelin (from MathWorks community)
%
%   Sample use: graphing standard deviations vs. means of climate output,
%   but colorscaled by latitude using 18 easily identifiable finite colors
%   in 180/18 = 10 degree intervals ([latyy] must be same size as [means], 
%   [stddev], and is usually a meshgridded version of the CMIP5 variable 
%   [lat]):
%   [lat_colors] = DATA_COLORMAP(latyy,18,'colormap','jet');
%   scatter(means,stddevs,10,lat_colors);
%
%   DATA_COLORMAP(...,'[flag]',[params],...) modify program run 
%   as below:
%       'integer'           - rounds colorbar split values down to the
%                             closest integer
%       'range',[rge]       - manually set range over which colormap should
%                             be set (i.e., if the data range is [-90,90],
%                             setting 'range',[-20,20] will set all values
%                             above 20 to the top colorbar color, and all
%                             values below -20 to the bottom). If [rge] =
%                             'def', nothing will change (default option
%                             setting kept to increase code flexibility)
%       'colormap',[cmap]   - set colormap manually. [cmap] can either be a
%                             string (valid MATLAB colormap name or name of
%                             a colormap in standard_colormaps.mat), a
%                             numeric RGB array, or a CSS3 supported color
%                             name . If manually inputted RGB array doesn't
%                             have a length divisble by [color_split],
%                             piece-wise linear interpolation is used to
%                             fill in colors. If fewer colors are provided
%                             than [color_split], only those colors are
%                             used and a warning is thrown.
%       'treat_nans',[log]  - whether to set NaNs to a specific color set
%                             below (by default true; not doing so can
%                             result in strange outlier points).
%       'nan_color',[num]   - set NaN color to which NaN members are
%                             assigned if treat_nans == true. By default
%                             [0.9 0.9 0.9]. Must be RGB value. 
%
%   Note: the file [standard_colormaps.mat], included with the code
%   package, includes a few useful colormaps for values that can take on
%   both negative and positive values, and can be used by setting the flag
%   'colormap' to one of the following names: 
%           - 'ColorRtB':   a colormap going from red to blue, passing 
%                           close to white at its center. 
%           - 'ColorRtWtB': a colormap passing from red to blue, but with a
%                           broad white center
%           - 'ColorRtGtB': a colormap passing from red to blue, but with a
%                           broad grey center (useful for scatterplots on a
%                           white background)
%
%   All directories listed as [____] are set in various_defaults.m
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last edit: 10/17/2017

%Vectorize dataset and make sure to avoid floating-point precision issues
dataset = round(dataset(:),8,'significant');

%Set defaults
cmap = flipud(jet(color_split));
clim_range = range(dataset);
start_value = min(dataset);
def_cmaps = matfile('standard_colormaps.mat');
treat_nans = true;
nan_color = [0.9 0.9 0.9];

%Set function modifiers
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'integer'}
                clim_range = ceil(max(dataset))-floor(min(dataset)); 
                start_value = floor(start_value);
            case {'range'}
                clim_range_tmp = varargin{in_idx+1}; varargin{in_idx+1} = 0;
                %Allow backwards compatibility for if range is given in
                %(...,#,#,...) format
                if length(clim_range_tmp) == 1; clim_range_tmp(2) = varargin{in_idx+2}; end
                %If 'range','def', do nothing (default option to increase
                %code flexibility)
                if ~isa(clim_range_tmp,'char') || ~strcmp(clim_range_tmp,'def')
                    clim_range = range(clim_range_tmp);
                    start_value = clim_range_tmp(1);
                end
                clear clim_range_tmp
            case {'colormap'}
                cmap = varargin{in_idx+1}; varargin{in_idx+1} = 0; %to prevent switch from tripping
                if isa(cmap,'char')
                    colormap_list = {'parula';'jet';'hsv';'hot';'cool';'spring';'summer';'autumn';'winter';'gray';'bone';'copper';'pink';'lines';'colorcube';'prism';'flag';'white'};
                    cust_colormap_list = fieldnames(def_cmaps);
                    if ~strcmp(cmap,'def')
                        if ~isempty(find(strcmp(cmap,cust_colormap_list),1))
                            cmap = def_cmaps.(cmap);
                        elseif isempty(find(strcmp(cmap,colormap_list),1))
                            %if a CSS3 color, can define using rgb
                            try cmap = rgb(cmap);
                            %otherwise throw error, since neither matlab or
                            %CSS3 supported colors
                            catch
                                error('DATA_COLORMAP:InvalidCmapName',['"',cmap,'" is not a valid MATLAB colormap or CSS3 color name; ',...
                                    'lists can be found <a href="http://www.mathworks.com/help/matlab/ref/colormap.html#input_argument_name">here</a> for MATLAB colormaps and ',...
                                    '<a href="https://www.w3.org/TR/css3-color/">here</a> for CSS3 colors.'])
                            end
                        else
                            cmap = eval([cmap,'(',num2str(color_split),')']); 
                        end
                    else
                        cmap = flipud(jet(color_split));
                    end
                end
                if size(cmap,1) < color_split
                    warning('DATA_COLORMAP:TooFewColors',['Inputted colormap has fewer colors than color_split, only ',num2str(size(cmap,1)),' color(s) will be used'])
                    color_split = size(cmap,1);
                end
            case {'treat_nans'}
                treat_nans = varargin{in_idx+1};
            case {'nan_color'}
                nan_color = varargin{in_idx+1}; varargin{in_idx+1} = 0;
        end
    end
end

%Modify colormap if it contains more colors than desired in [color_split]
if size(cmap,1) ~= color_split
    if mod(size(cmap,1),color_split)~=0 %If interpolation needed (# of colors not an integer multiple of [color_split])
        %Use piecewise linear interpolation to create new colormap with the
        %desired number of colors
        colorscale_split = zeros(color_split,3);
        for rgb_idx = 1:3
            colorscale_split(:,rgb_idx) = piecelin((1:size(cmap,1))',cmap(:,rgb_idx),(1:size(cmap,1)/(color_split):size(cmap,1))');
        end
        cmap = colorscale_split; clear colorscale_split
    else %If number of colors is just integer multiple of [color_split]
        cmap = cmap(1:size(cmap,1)/(color_split):size(cmap,1),:);
    end
end


%Set #colors data values to round data to
color_steps = start_value+(0:color_split-1)*(clim_range/(color_split));

%Assign data to colors from cmap by finding the closest color_steps value
%to each data point by minimizing distance to the midpoint of the 'color
%step' defined as [color_steps(k),color_steps(k+1)).
rgb_adj = zeros(length(dataset),3);
for i = 1:length(dataset)
    if size(cmap,1) == 1 %If just one color inputted, assign that to all data values
        rgb_adj(i,:) = cmap;
    else
        if treat_nans && isnan(dataset(i))
            rgb_adj(i,:) = nan_color;
        else
        %if ~isnan(dataset(i)) && treat_nans
            [~,min_idx] = min(abs(dataset(i)-(color_steps+0.5*abs(color_steps(2)-color_steps(1)))));
            rgb_adj(i,:) = cmap(min_idx,:);
        %else
            %rgb_adj(i,:) = nan_color;
        end
    end
end

%Output data
varargout{1} = rgb_adj;
if nargout >= 2
    varargout{2} = cmap;
    if nargout == 3
        varargout{3} = [start_value start_value+clim_range];
    end
end
end