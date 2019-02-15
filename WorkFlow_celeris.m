simStart = '2019-01-14T13.00.00.Z'; % hindcast date in Z
simEnd = '2019-01-14T14.00.00.Z'; % end date for batchruns 
homeDir = "C:\Users\Waves.DESKTOP-DRJMBVD\Desktop\Celeris_dump\Celeris_Duck\Celeris_Duck";
sim_time=4*60; % wall clock time (in seconds)  should be 
hourDT = 1; % run model every hour between start and end


%% run 
% generate list of forecast_dates
endSimdatetime = datetime(simEnd, 'InputFormat','yyyy-mm-dd''T''HH.mm.ss''.Z');
startSimdatetime = datetime(simStart, 'InputFormat','yyyy-mm-dd''T''HH.mm.ss''.Z');
simList = startSimdatetime:hours(hourDT):endSimdatetime;
for ii=1:length(simList)
   simulation_date = datestr(simList(ii), 'yyyy-mm-ddTHH.MM.ss.Z');
    cd(homeDir)
    DuckSurf_windwave_nowcast(simulation_date, sim_time)    
end