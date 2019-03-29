function [WL] = getWLFRF(d1, d2, gnum)
% getwL retrieves Water Level data from the thredds server 
%   INPUTS 
%       d1-start date in matlab datenum format - ex. datenum(2015,10,2)
%       d2-end date in matlab datenum format - see above
%       gnum is the gauge number (11 only for WL)
%           1 = gauge 11 (end of pier) - all that is currently available
%
%altered by Julia Fiedler, Jan 20, 2019
%jfiedler@ucsd.edu
%MATLAB v.2018b
%TODO: Check input on time, if datenum or datetime format and proceed
%accordingly. Optimal: datetime, probably.
%
% % % % % % % % 
%% Main code 
% defining data location url (1st part)
svrloc='https://134.164.129.55/thredds/dodsC/FRF';  % The prefix for the CHL thredds server
svrloc='https://chldata.erdc.dren.mil/thredds/dodsC/frf';  % The prefix for the CHL thredds server

% defining 2nd part
%TODO: Can change which gauge to point to here

if gnum==1
    urlback='/oceanography/waterlevel/eopNoaaTide'; % Water Level 
elseif gnum==2
    urlback='/oceanography/waterlevel/eopNoaaTide';    
else
    disp 'please visit http://chldata.erdc.dren.mil/thredds/catalog/frf/catalog.html\n' ...
          'and browse to the gauge of interest and select proper url and add to program'
end

d1vec = datevec(d1);
d2vec = datevec(d2);
%% set URL 
if d1vec(1:2) == d2vec(1:2) %if the year + month are the same, use the url for only that month
    urlyear = num2str(d1vec(1));
    urlmonth = num2str(d1vec(2),'%02.f');
    url = [svrloc urlback '/' urlyear '/FRF-ocean_waterlevel_eopNoaaTide_' urlyear urlmonth '.nc'];
else
	endUrl = '/eopNoaaTide.ncml';  % if not use the generic url
	url=strcat(svrloc,urlback,endUrl); % combining first and 2nd part of url string
end
    
%% go now find index  
% pulling down time now
time=ncread(url,'time'); % downloading time from server
tunit= ncreadatt(url,'time','units'); % reading attributes of variable time
% converting time to matlab datenum
mtime=time/(3600.0*24)+datenum(1970,1,1);
% finding index that corresponds to dates of interest
sprintf('WL record starts %s and ends %s', datestr(min(mtime)),datestr(max(mtime)));
itime=find(d1 < mtime & d2> mtime); % indicies in netCDF record of data of interest
%% pull data with appropriate time index 
% pulling data from server with itime index 
WL.time = mtime(itime);        % record of time in matlab datetime 
WL.WL = ncread(url,'waterLevel',min(itime),length(itime));  % measured water level
WL.PredictedWL = ncread(url,'predictedWaterLevel',min(itime),length(itime));  % predicted WL
WL.name = ncread(url,'station_name'); % station name 
WL.lat = ncread(url,'latitude'); % latitutde 
WL.lon = ncread(url, 'longitude'); % lon 

end