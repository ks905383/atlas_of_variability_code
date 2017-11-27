function varargout = subplot_pos(num_subs)
% SUBPLOT_POS   get standardized subplot positions.
%
%   position_set = SUBPLOT_POS(num_subs) outputs a num_subs x 4 array of
%   position values [position_set] for a 2 x mod([num_subs]/2) grid of
%   subplots to use for subplot positioning.
%
%   [position_set,pointsx,pointsy] = SUBPLOT_POS(num_subs) also outputs the
%   width and height of image in points (1/72th of an inch) [pointsx] and
%   [pointsy] for which the subplot position values were optimized.
%
%   Currently only supports up to 6 subplots. 
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last edit: 09/16/2016

if num_subs ==1 
    position_set(1,1:4) = [0.05,0.05,0.85,0.90];
    pointsy = 576; pointsx = 792;
    colorbar_pos = [];
    position_set_alt = position_set;
elseif num_subs == 2
    position_set(1,1:4) = [0.235,0.55,0.45,0.4];
    position_set(2,1:4) = [0.235,0.10,0.45,0.4];
    pointsy = 576; pointsx = 792;
    colorbar_pos = [];
    position_set_alt = position_set;
elseif num_subs <= 4 && num_subs > 2
    position_set(1,1:4) = [0.01,0.55,0.45,0.4];
    position_set(2,1:4) = [0.46,0.55,0.45,0.4];
    position_set(3,1:4) = [0.01,0.10,0.45,0.4];
    position_set(4,1:4) = [0.46,0.10,0.45,0.4];
    pointsy = 576; pointsx = 792;
    colorbar_pos = position_set+[0,0.013,-0.4235,-0.035;0.425,0.013,-0.4235,-0.035;0,0.013,-0.4235,-0.035;0.425,0.013,-0.4235,-0.035];
    position_set_alt = position_set;
    position_set_alt([1 3],:) = position_set_alt([1 3],:) + [0.05 0 -0.01 0; 0.05 0 -0.01 0];
    position_set_alt([2 4],:) = position_set_alt([2 4],:) + [0.04 0 -0.01 0; 0.04 0 -0.01 0];
elseif num_subs <= 6 && num_subs > 4
    position_set(1,1:4) = [0.05,0.6667,0.43,0.25];
    position_set(2,1:4) = [0.48,0.6667,0.43,0.25];
    position_set(3,1:4) = [0.05,0.3833,0.43,0.25];
    position_set(4,1:4) = [0.48,0.3833,0.43,0.25];
    position_set(5,1:4) = [0.05,0.10,0.43,0.25];
    position_set(6,1:4) = [0.48,0.10,0.43,0.25];
    pointsy = 792; pointsx = 576;
    colorbar_pos = [];
    position_set_alt = position_set;
elseif num_subs > 6
    error('subplot_pos:NotSupported',['subplot_pos currently only supports up to 6 subplots, ',num2str(num_subs),' (inputted) is too many.'])
end

varargout{1} = position_set;
if nargout > 1
    varargout{2} = pointsx;
    varargout{3} = pointsy;
    varargout{4} = colorbar_pos;
    varargout{5} = position_set_alt;
end

end