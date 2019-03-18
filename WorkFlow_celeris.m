%% Celeris master Work flow for Coastal Model Test Bed
%    This script can operate in hindcast mode by inputing start and end 
%    dates or in Nowcast mode modeling the most recent day
%    
%   INPUTS:
%       simStart - time to start a batch of simulation(s).  If left blank (ie. []) 
%           then will model the most recent day.  To model in hindcast input
%           time in format of '2019-01-24T13.00.00.Z'.  All times are in
%           UTC to match server time. 
%       simEnd - time to end a batch of simulation(s).  If left blank (ie. []) 
%           then will model the most recent day.  To model in hindcast input
%           time in format of '2019-01-24T13.00.00.Z'.  All times are in
%           UTC to match server time. 
%       homeDir - this is where the celeris model exe and all
%           supplimental files are located.  must also have the yaml_files
%           directory to make netCDF files
%       hourDT - this is the increment to run the model at (default is 1,
%           for every hour)
%
%
%
%% set inputs 
simStart = [];   % hindcast date in Z
simEnd =   [];   % end date for batchruns 
% set home working directory 
homeDir = "D:\CMTB_Celeris";
% wall clock time (in seconds)  should be 30, will save the last 17 minutes
sim_time=30*60; 
% run model every hour between start and end  (how often to run between dates)
hourDT = 1; 
%% run model loop
% generate list of forecast_dates
if isempty(simStart);
    endSimdatetime = dateshift(datetime('now'), 'start', 'hour');
    startSimdatetime = endSimdatetime - hours(24);
    hourDT = 1;
else;
    endSimdatetime = datetime(simEnd, 'InputFormat','yyyy-MM-dd''T''HH.mm.ss''.Z');
    startSimdatetime = datetime(simStart, 'InputFormat','yyyy-MM-dd''T''HH.mm.ss''.Z');
    hourDT = hourDT; 
end
% generate simulation list 
simList = startSimdatetime:hours(hourDT):endSimdatetime;
for ii=1:length(simList)
    simulation_date = datestr(simList(ii), 'yyyy-mm-ddTHH.MM.ss.Z');
    cd(homeDir)
    DuckSurf_windwave_nowcast(simulation_date, sim_time)    
end