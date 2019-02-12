function load_FRFwave_lidar(fname)

waveTime=ncread(fname,'time');
waveHs=ncread(fname,'waveHs');
waveTp=ncread(fname,'waveTp');
waveTm=ncread(fname,'waveTm');
xFRF=ncread(fname,'xFRF');
yFRF=ncread(fname,'yFRF');
waterLevel=ncread(fname,'waterLevel');
waveFrequency=ncread(fname,'waveFrequency');
waveHsIG=ncread(fname,'waveHsIG');
waveHsTotal=ncread(fname,'waveHsTotal');
waveAsymmetry=ncread(fname,'waveAsymmetry');

eval(['save FRFwave_forecast_lidar.mat waveTime waveHs waveTp waveTm xFRF yFRF waterLevel waveFrequency waveHsIG waveHsTotal waveAsymmetry'])
