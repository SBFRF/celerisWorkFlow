%% Celeris master Work flow for Coastal Model Test Bed
%    This script can operate in hindcast mode by inputing start and end 
%    dates or with 

simStart = '2019-01-24T13.00.00.Z'; % hindcast date in Z
simEnd = '2019-01-24T14.00.00.Z'; % end date for batchruns 
% set home working directory 
homeDir = "C:\Users\Waves.DESKTOP-DRJMBVD\Desktop\Celeris_dump\Celeris_Duck\Celeris_Duck";
% wall clock time (in seconds)  should be 30, will save the last 17 minutes
sim_time=30*60; 
% run model every hour between start and end  (how often to run between dates)
hourDT = 1; 
%% run 
% generate list of forecast_dates
if isempty(simStart);
    endSimdatetime = dateshift(datetime('now'), 'start', 'hour');
    startSimdatetime = endSimdatetime - hours(24);
else;
    endSimdatetime = datetime(simEnd, 'InputFormat','yyyy-MM-dd''T''HH.mm.ss''.Z');
    startSimdatetime = datetime(simStart, 'InputFormat','yyyy-MM-dd''T''HH.mm.ss''.Z');
end
% generate simulation list 
simList = startSimdatetime:hours(hourDT):endSimdatetime;
for ii=1:length(simList)
    simulation_date = datestr(simList(ii), 'yyyy-mm-ddTHH.MM.ss.Z');
    cd(homeDir)
    DuckSurf_windwave_nowcast(simulation_date, sim_time)    
end