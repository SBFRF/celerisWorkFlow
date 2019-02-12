% loads to-file output from Celeris

load time_axis.txt -ascii
load array.txt -ascii

x_ind=array(:,1);
y_ind=array(:,2);
etacol=array(:,3);

numx=max(x_ind)-min(x_ind)+1;
numy=max(y_ind)-min(y_ind)+1;

min_x_inst=min(instruments(:,1));
max_x_inst=max(instruments(:,1));
min_y_inst=min(instruments(:,2));
max_y_inst=max(instruments(:,2));
dx_inst=(max_x_inst-min_x_inst)/(numx-1);
x_inst=[min_x_inst:dx_inst:max_x_inst]+xoffset;
y_inst=instruments(1,2)+yoffset;

nt=length(time_axis);
eta=zeros(numx,numy,nt);
block=numx*numy;
for n=1:nt
    bs=1+(n-1)*block;
    be=n*block;
    eta(:,:,n)=reshape(etacol(bs:be),[numx,numy]);
end

 for n=1:2710 %nt
     plot(x_inst,squeeze(eta(:,1,n))'-ho')
     axis([-100 110 -Inf Inf])
     pause(.1)
 end
 ho=squeeze(eta(:,1,1));
 Total_depth=squeeze(eta(:,1,:))'-ho';
 shore_ie=find(ho<=0.001,1);
 min_depth=0.01;
 for n=1:2710
     for i=shore_ie:-1:1
         if Total_depth(i)<min_depth
             runup(n)=squeeze(eta(i,1,n))';
             break
         end
     end
 end
 
 plot(runup)
 
 
 

Hs=zeros(numx,numy);
zmean=Hs;
Tp=Hs;
Tm=Hs;

nt_start=round(nt/6);
for i=1:numx
    for j=1:numy
        [Hs(i,j),zmean(i,j),Tp(i,j),Tm(i,j),f,Se,f_ave,S_ave]=spectral(time_axis(nt_start:nt),squeeze(eta(i,j,nt_start:nt))');
    end
end


% 8m-array processing
inst_ind=1;
x_inst_FRF(inst_ind)=917;
[Hmo_cut_spectrum(inst_ind),Hmo_all(inst_ind),Tp_all(inst_ind)]=load_FRFinst(fname_writeout8,1,forecast_num_FRF);

% AWAC 45 processing
inst_ind=2;
x_inst_FRF(inst_ind)=400;
[Hmo_cut_spectrum(inst_ind),Hmo_all(inst_ind),Tp_all(inst_ind)]=load_FRFinst(fname_writeout45,1,forecast_num_FRF);

% % ADOP 35 processing
% inst_ind=3;
% x_inst_FRF(inst_ind)=300;
% [Hmo_cut_spectrum(inst_ind),Hmo_all(inst_ind),Tp_all(inst_ind)]=load_FRFinst(fname_writeout35,1);
% 
% % XP 125m processing
% inst_ind=4;
% x_inst_FRF(inst_ind)=125;
% [Hmo_cut_spectrum(inst_ind),Hmo_all(inst_ind),Tp_all(inst_ind)]=load_FRFinst(fname_writeout19,2);

% Lidar Hydro
fname_lidar='FRF-ocean_waves_lidarHydrodynamics_201901.nc'
load_FRFwave_lidar(fname_lidar);
load FRFwave_forecast_lidar.mat
Nt_lidar=find(waveTime>=forecast_num_FRF,1);
time_reference = datenum('1970', 'yyyy');
time_EDT = time_reference + double(waveTime(Nt_lidar))/24/60/60-5/24;  % EDT time
str1=['Nowcast Time: ' datestr(time_EDT,'yyyy-mm-dd HH:MM') ' EDT'];

% Lidar Runup
fname_lidar='FRF-ocean_waves_lidarRunup_201901.nc'
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











