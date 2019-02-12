function load_FRFwave(fname)

waveTime=ncread(fname,'time');
Hs=ncread(fname,'waveHs');
Tp=ncread(fname,'waveTp');
Dp=ncread(fname,'waveMeanDirection');
nominaldepth=ncread(fname,'nominalDepth');
wavedepth=ncread(fname,'depth');
waveFrequency=ncread(fname,'waveFrequency');
waveMeanDirection=ncread(fname,'waveDirectionBins');
waveEnergyDensity=ncread(fname,'directionalWaveEnergyDensity');

eval(['save FRFwave_forecast.mat waveTime Hs Tp Dp nominaldepth wavedepth waveFrequency waveMeanDirection waveEnergyDensity'])
