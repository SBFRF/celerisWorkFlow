% loads to-file output from Celeris
% written by Pat Lynett, USC and modified Spicer Bak, USACE 

load repeat_time.txt -ascii  % repeat_time in minutes
load time_axis.txt -ascii
load array.txt -ascii
freqInterp = 0.5;       % what resolution do we need to sample the timeseries output at 
tEnd =17.0667;          %[m] always Take the last 17 minutes of runtime (determined by repeat time)

% parse out spatial time series
x_ind=array(:,1);
y_ind=array(:,2);
etacol=array(:,3);
pcol = array(:, 4); % total water depth * u 
qcol = array(:, 5); % total water depth * v
timeArray = time_axis(:); % parse out time

%% parse into dimensioned arrays
numx=max(x_ind)-min(x_ind)+1;
numy=max(y_ind)-min(y_ind)+1;
nt=length(time_axis);
eta=zeros(numx,numy,nt);
p = zeros(numx,numy,nt);
q = zeros(numx,numy,nt);
block=numx*numy;
for n=1:nt
    bs=1+(n-1)*block;
    be=n*block;
    eta(:,:,n)=reshape(etacol(bs:be),[numx,numy]);
    p(:,:,n)=reshape(pcol(bs:be), [numx, numy]);
    q(:,:,n)=reshape(qcol(bs:be), [numx, numy]);
end

%% load instrument positions
min_x_inst=min(instruments(:,1));
max_x_inst=max(instruments(:,1));
min_y_inst=min(instruments(:,2));
max_y_inst=max(instruments(:,2));
dx_inst=(max_x_inst-min_x_inst)/(numx-1);
x_inst=[min_x_inst:dx_inst:max_x_inst]+xoffset;
y_inst=instruments(1,2)+yoffset;
%% calculate runup from model output 
ho=squeeze(eta(:,1,1));  % original Depth 
Total_depth=squeeze(eta(:,1,:))'-ho';
shore_ie=find(ho<=0.001,1);
min_depth=0.01;
runup=zeros(nt,1);
for n=1:nt
    for i=shore_ie:-1:1
        if Total_depth(i)<min_depth
            runup(n)=squeeze(eta(i,1,n))';
            break
        end
    end
end
runup_trunc=runup(:) + water_level_change;
% this function seems to bomb at times, we will save as NaN, but save time
% series to calculate later
try
    [R2]=calc_runup_dist(runup_trunc);
catch
    [R2]=NaN;
end

ho_dry=ho(1:shore_ie) + water_level_change;
xFRF_dry=x_inst(1:shore_ie);

Hs=zeros(numx,numy);
zmean=Hs;
Tp=Hs;
Tm=Hs;

nt_start=find(time_axis>time_axis(nt)-repeat_time*60,1);
for i=1:numx
    for j=1:numy
%       [Hs(i,j),zmean(i,j),Tp(i,j),Tm(i,j),f,Se,f_ave,S_ave]=spectral(outTime(idxStart:end), squeeze(eta(i,j,idxStart:end))');
        [Hs(i,j),zmean(i,j),Tp(i,j),Tm(i,j),f,Se,f_ave,S_ave]=spectral(time_axis(nt_start:nt),squeeze(eta(i,j,nt_start:nt))');
    end
end
%% now Interpolate time series to regular output 
%     This forces data to be the same length for netCDF files 

outTime = 1:freqInterp:timeArray(end);
if outTime(end)/60 > tEnd
    idxStart = length(outTime)-ceil(tEnd*60/freqInterp)+1;
else % if the simulation is shorter than tEnd
    idxStart = 1; 
end 
etaInterp = interp2(squeeze(x_inst), squeeze(timeArray), squeeze(eta)',  x_inst, outTime', 'spline');
uInterp = interp2(squeeze(x_inst), squeeze(timeArray), squeeze(p)', x_inst, outTime', 'spline'); % total water depth * u 
vInterp = interp2(squeeze(x_inst), squeeze(timeArray), squeeze(q)', x_inst, outTime', 'spline'); % total water depth * v
runupInterp = interp1(squeeze(timeArray), squeeze(runup), outTime', 'spline');
% if model grid resolution changes, we may need to interpolate in x/y 
%% Write NetCDF output

if exist('netCDFcode', 'dir')
    addpath(genpath('netCDFcode'))
else
    error('netCDF package not Found!  Please add netCDFcode to your search path');
end
if exist('yaml_files', 'dir')
    addpath(genpath('yaml_files'))
else
    error('netCDF package not Found!  Please add netCDFcode to your search path');
end

% create data structure for netCDF file output
dataIn.time = epoch2Matlab(waveTime(Nt));   % this should be matlab time 
dataIn.tsTime = (squeeze(outTime(idxStart:end)) - outTime(idxStart))';
dataIn.eta = (water_level_change + etaInterp(idxStart:end, : ))'; % permute(water_level_change + etaInterp(idxStart:end, : ), [3, 1, 2]);  
dataIn.velocityU = (uInterp(idxStart:end, :))';
dataIn.velocityV = (vInterp(idxStart:end, :))';
dataIn.xFRF = x_inst;
dataIn.yFRF = y_inst;
dataIn.station_name = "celeris Model Profile";
dataIn.totalWaterLevel = permute(R2, [2,1]) ;
dataIn.totalWaterLevelTS =  runupInterp(idxStart:end); % time series
% note check for what makes dimensions, why am i getting weird sizes 
matlab2netCDF(dataIn, globalYamlFileName, varYamlFileName, 1, fnameOut);
%% process for gauge comparisons 
addTime  =  3600;  % time in seconds to determine data gathering window
dataCollectStart = datenum(datetime(waveTime(Nt), 'convertfrom','posixtime'));
dataCollectionEnd = datenum(datetime(waveTime(Nt)+addTime, 'convertfrom','posixtime')); % add 1 hour 

% loop over inst_ind  associated with getwaveFRF 
% index = 1:6;  %indicies of instruments to compare
try
    % 8m-array processing
    % x_inst_FRF(inst_ind)=917;
    inst_ind=1;
    array8m = getwaveFRF(dataCollectStart, dataCollectionEnd, 10);   % this could even get looped over
    [~, ~, ~, ~, yFRF, x_inst_FRF(inst_ind)] = frfCoord(array8m.lat, array8m.lon);  % get FRF coordinates for the gauge
    [Hmo_cut_spectrum(inst_ind),Hmo_all(inst_ind),Tp_all(inst_ind)]=load_FRFinst(array8m, waveTime(Nt));
catch
end
try
    % AWAC 45 processing
    % x_inst_FRF(inst_ind)=400;
    inst_ind=2;
    awac45m = getwaveFRF(dataCollectStart, dataCollectionEnd, 8);   % this could even get looped over
    [~, ~, ~, ~, yFRF, x_inst_FRF(inst_ind)] = frfCoord(awac45m.lat, awac45m.lon);  % get FRF coordinates for the gauge
    [Hmo_cut_spectrum(inst_ind),Hmo_all(inst_ind),Tp_all(inst_ind)]=load_FRFinst(awac45m, waveTime(Nt));
catch
end
try
    % % ADOP 35 processing
    % x_inst_FRF(inst_ind)=300;
    inst_ind=3;
    adop35m = getwaveFRF(dataCollectStart, dataCollectionEnd, 8);   % this could even get looped over
    [~, ~, ~, ~, yFRF, x_inst_FRF(inst_ind)] = frfCoord(adop35m.lat, adop35m.lon);  % get FRF coordinates for the gauge
    [Hmo_cut_spectrum(inst_ind),Hmo_all(inst_ind),Tp_all(inst_ind)]=load_FRFinst(adop35m,waveTime(Nt));
catch
end
%
% % XP 125m processing
% inst_ind=4;
% x_inst_FRF(inst_ind)=125;
% [Hmo_cut_spectrum(inst_ind),Hmo_all(inst_ind),Tp_all(inst_ind)]=load_FRFinst(fname_writeout19,2);
try
    % Lidar Hydro
    fname_lidar='FRF-ocean_waves_lidarHydrodynamics_201901.nc';
    load_FRFwave_lidar(fname_lidar);
    load FRFwave_forecast_lidar.mat
    Nt_lidar=find(waveTime>=forecast_num_FRF,1);
    time_reference = datenum('1970', 'yyyy');
    time_EDT = time_reference + double(waveTime(Nt_lidar))/24/60/60-5/24;  % EDT time
    str1=['Nowcast Time: ' datestr(time_EDT,'yyyy-mm-dd HH:MM') ' EDT'];
    
    % Lidar Runup
    fname_lidar='FRF-ocean_waves_lidarRunup_201901.nc';
    load_FRFwave_lidar(fname_lidar);
    load FRFwave_forecast_lidar.mat
    
    hf6=figure(6);
    clf
    plot(x_inst,Hs(:,1), x_inst,1*(zmean(:,1)+water_level_change))
    hold on
    plot(x_inst_FRF,Hmo_cut_spectrum,'ro')
    plot(x_inst_FRF,Hmo_all,'go')
    plot(xFRF,waveHsTotal(:,Nt_lidar),'g-','LineWidth',2)
    plot(xFRF,1*waterLevel(:,Nt_lidar),'b-','LineWidth',2)
    title('Model-Data Comparison along y_{FRF}=940 m')
    xlabel('x_{FRF} [m]')
    ylabel('Elevation [m]')
    legend('Modeled H_s','Modeled setup x 10','Observed H_s (onshore only, >6s)','Observed H_s (all)','Lidar Hs','Lidar Setup','Location','SouthEast')
catch
end










