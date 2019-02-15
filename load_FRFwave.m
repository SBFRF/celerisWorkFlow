function load_FRFwave(fname)
try  % loads from flat file
waveTime=ncread(fname,'time');
Hs=ncread(fname,'waveHs');
Tp=ncread(fname,'waveTp');
Dp=ncread(fname,'waveMeanDirection');
nominaldepth=ncread(fname,'nominalDepth');
wavedepth=ncread(fname,'depth');
waveFrequency=ncread(fname,'waveFrequency');
waveMeanDirection=ncread(fname,'waveDirectionBins');
waveEnergyDensity=ncread(fname,'directionalWaveEnergyDensity');
catch  % this will parse a structure from the getWaveFRF, from thredds 
    waveTime = fname.time;
    Hs = fname.Hs;
    Tp = fname.Tp;
    Dp = fname.Dp;
    nominaldepth = median(fname.depth);  % incase there are multiple values
    wavedepth = fname.depth;
    waveFrequency = 1/fname.Tp;
    waveMeanDirection = fname.dirpeak;  % this isn't totally right
    waveEnergyDensity = fname.spec2D;
end
eval(['save FRFwave_forecast.mat waveTime Hs Tp Dp nominaldepth wavedepth waveFrequency waveMeanDirection waveEnergyDensity'])
