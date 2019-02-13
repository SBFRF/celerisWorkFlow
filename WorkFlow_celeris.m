forecast_date= '2019-01-14T13.00.00.Z'; % hindcast date in Z
homeDir = "C:\Users\Waves.DESKTOP-DRJMBVD\Desktop\Celeris_dump\Celeris_Duck\Celeris_Duck";
sim_time=5*60; % wall clock time (in seconds) 


%% run 
cd(homeDir)
DuckSurf_windwave_nowcast(forecast_date, sim_time)