function [lat_idxs,lon_idxs] = subset_find(lat,lon,lat_alt,lon_alt)
% SUBSET_FIND   get indices of a block of lat/lon values within a larger
%               grid.
%
%   [lat_idxs,lon_idxs] = SUBSET_FIND(lat,lon,lat_alt,lon_alt) returns the
%   (subscript) indices of [lat] and [lon] corresponding to the subset of
%   the lat x lon grid represented by [lat_alt] and [lon_alt]. The inputs
%   can either be vectors (nlat x 1, nlon x 1, etc.) or arrays (lat = [nlat
%   x nlon], etc.). 
%
%   Example: given a 2-D grid with (lat, lon) coordinates (a:d, 1:8):
%           a1 a2 a3 a4 a5 a6 a7 a8
%           b1 b2 b3 b4 b5 b6 b7 b8
%           c1 c2 c3 c4 c5 c6 c7 c8
%           d1 d2 d3 d4 d5 d6 d7 d8
%       if inputting a (lat_alt, lon_alt) subset of this grid with values:
%                 b3 b4
%                 c3 c4
%                 d3 d4
%       SUBSET_FIND(lat,lon,lat_alt,lon_alt) would return 
%           lat_idxs = [2 3 4];
%           lon_idxs = [3 4];
%
%   Sample uses: if (for, say a climate model) a nlon x nlat array of area
%   weights exists for a global grid, but some quantity is saved on only a
%   subset of that grid (say, North America), subset_find allows for a
%   quick subsetting of the area weight array to match the corresponding
%   subsetted grid, as long as a vector (or 2d array) of lat lon values
%   exists for the area weight grid and lat_alt lon_alt for the subset
%   grid. 
%
%   More generally, SUBSET_FIND allows the indices of any sub-matrix to be
%   found in any matrix. 
%
%   See also 
%
%   For questions/comments, contact Kevin Schwarzwald
%   kschwarzwald@uchicago.edu
%   Last modified: 01/25/2017


%lat/lon = full nlon x nlat grid, to be shortened to match
%lat_alt/lon_alt = clipped, corresponding to the smaller one

if ~isequal([size(lat) size(lon)],[size(lat_alt) size(lon_alt)])
    %Find start of overlap, if same grid
    test_idx_lat = find(lat == lat_alt(1,1));
    if ~isempty(test_idx_lat)
        %For all equals, try to find the one that's actually at
        %the start of the overlapping block. 10/9 are
        %arbitrary as length of test overlap
        for i = 1:length(test_idx_lat)
            if length(lat_alt)>10
                if length(lat)-max(test_idx_lat)>9 || (size(lat,2)>1 && numel(lat)-max(test_idx_lat)>9)
                    test_idx_max = 9;
                else %If cannot get 10 consecutive pixels to test, get largest subset
                    test_idx_max = length(lat)-max(test_idx_lat)-1;
                    warning('SUBSET_FIND:SmallTest',['Tested for convergence over ',num2str(test_idx_max),' latitudes, hopefully that"s enough.'])
                end
            else
                test_idx_max = length(lat_alt)-1;
            end
            %Find if this index anchors a block of
            %overlapping lat/lon grid values
            if isequal(reshape(lat(test_idx_lat(i):test_idx_lat(i)+test_idx_max),[test_idx_max+1,1]),...
                    reshape(lat_alt(1:test_idx_max+1,1),[test_idx_max+1,1]))
                start_idx_lat = i;
            end
        end
    else
        error('SUBSET_FIND:NoLatMatch',['The latitude ',num2str(lat_alt(1,1)),...
            ' is not present in the full latitude array, no match/subset possible.'])
    end
    %If vector lat/lon grid, repeat analysis to find lon start
    %index as well
    if size(lat_alt,2) == 1
        test_idx_lon = find(lon == lon_alt(1,1));
        clear start_idx
        if ~isempty(test_idx_lon)
            %For all equals, try to find the one that's actually at
            %the start of the overlapping block. 10/9 are
            %arbitrary as length of test overlap
            for i = 1:length(test_idx_lon)
                if length(lon_alt)>10
                    if length(lon)-max(test_idx_lon)>9
                        test_idx_max = 9;
                    else %If cannot get 10 consecutive pixels to test, get largest subset
                        test_idx_max = length(lon)-max(test_idx_lon)-1;
                        warning('SUBSET_FIND:SmallTest',['Tested for convergence over ',num2str(test_idx_max),' longitudes, hopefully that"s enough.'])
                    end
                else
                    test_idx_max = length(lon_alt)-1;
                end
                %Find if this index anchors a block of
                %overlapping lat/lon grid values
                if isequal(reshape(lon(test_idx_lon(i):test_idx_lon(i)+test_idx_max),[test_idx_max+1,1]),...
                        reshape(lon_alt(1:test_idx_max+1,1),[test_idx_max+1,1]))
                    start_idx_lon = i;
                end
            end
        else
            error('SUBSET_FIND:NoLonMatch',['The longitude ',num2str(lon_alt(1,1)),...
                ' is not present in the full longitude array, no match/subset possible.'])
        end
    end
    %If this overlap exists, use it as a seed to get new indices to clip
    if exist('start_idx_lat','var')
        if size(lat,2) > 1
            [strt_x,strt_y] = ind2sub(size(lat),test_idx_lat(start_idx_lat));
            lat_idxs = strt_y:strt_y+size(lat_alt,2)-1;
            lon_idxs = strt_x:strt_x+size(lat_alt,1)-1;
        else
            lat_idxs = test_idx_lat(start_idx_lat):test_idx_lat(start_idx_lat)+size(lat_alt,1)-1;
            lon_idxs = test_idx_lon(start_idx_lon):test_idx_lon(start_idx_lon)+size(lon_alt,1)-1;
        end
    else
        error('SUBSET_FIND:NoGridMatch',['No block of >= ',num2str(test_idx_max),...
            ' matching latitude/longitude values were found between the full and alternate lat/lon grid inputs.'])
    end
else
    warning('SUBSET_FIND:SameGrid','The full and alternate lat/lon grids inputted are identical. No subset made.')
    if size(lat,2) > 1
        lat_idxs = 1:size(lat,2);
        lon_idxs = 1:size(lat,1);
    else
        lat_idxs = 1:length(lat);
        lon_idxs = 1:length(lon);
    end
end