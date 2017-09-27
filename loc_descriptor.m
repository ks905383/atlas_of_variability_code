% LOC_DESCRIPTOR    reverse geocoding from lat,lon pair using Google Maps
%                   API tools
%   loc_desc = LOC_DESCRIPTOR(lat,lon) returns a string giving the most
%   complete form of the location's address stored in the Google Maps
%   database.
%
%   Description of API can be found <a href =
%   "https://developers.google.com/maps/documentation/geocoding/intro?hl=de">here</a>
%
%   NOT YET IMPLEMENTED, BUT TO IMPLEMENT:
%       - set which level of address to spit out
%       - allow location ID as input? Maybe?
%
%   LOC_DESCRIPTOR(...,'flags',[params],...) modify program run as follows:
%       'api_key',[char]        - input api key (needed for more complex
%                                 reverse geocoding)
%       'result_type',[char]    - pipes into "result_type" parameter, i.e.
%                                 'political' for only political divisions
%                                 (country, state, county, etc.) (requires
%                                 api_key)
%
%   Last modified 06/21/2016 by Kevin Schwarzwald
%   For questions, comments, contact kschwarzwald@uchicago.edu


function loc_desc = loc_descriptor(lat,lon,varargin)
%% Set defaults and optional function flag behavior
api_key = [];
result_type = [];
cust_input = [];
%Set behavior of optional function flags
if (~isempty(varargin))
    for in_idx = 1:length(varargin)
        switch varargin{in_idx}
            case {'api_key'}
                api_key = varargin{in_idx+1};
            case {'result_type'}
                result_type = varargin{in_idx+1};
        end
    end
end

%% Read API Url
%Make sure lon is in [-180,180] format
if lon > 180; lon = -180 + (lon-180); end
%Make sure lat is in [-90,90] format
if lat > 90 || lat < -90; error('LOC_DESCRIPTOR:BadLat',['Latitude is inputted as ',...
        num2str(lat),', must be between -90 and 90.'])
end
%Get base API url string
googleapi_url = ['https://maps.googleapis.com/maps/api/geocode/json?latlng=',num2str(lat),',',num2str(lon)];

%Add result type if desired (i.e. 'political')
if ~isempty(result_type);googleapi_url = [googleapi_url,'&result_type=',result_type]; end

%Add any custom extension to the url
googleapi_url = [googleapi_url,cust_input];

%Add API key if provided (needed for more complex requests)
if ~isempty(api_key);googleapi_url = [googleapi_url,'&key=',api_key];end

%Get API output text from google maps
api_text = urlread(googleapi_url);

%% Various Error Messaegs
if strfind(api_text,'ZERO_RESULTS')
    error('LOC_DESCRIPTOR:ZeroResults',...
        ['Zero results in Maps API for the search ',googleapi_url(37:end),...
        ' (a possible reason could be that the desired point is in international waters)']);
end
if strfind(api_text,'REQUEST_DENIED')
   error('LOC_DESCRIPTOR:RequestDenied',...
       ['Error message shown: ',api_text(strfind(api_text,'"error_message"')+18:strfind(api_text,',')-1)])
end

%% Output
%Find easily readable address version of the API output text
form_add_idx = strfind(api_text,'"formatted_address"');
geo_idx = strfind(api_text,'"geometry"');

%Output readable text
loc_desc = api_text(form_add_idx(1)+23:geo_idx-1);

%Get rid of quotes, trailing comma
aux_find = regexp(loc_desc,'\,');
loc_desc = loc_desc(1:aux_find(end)-2);

end