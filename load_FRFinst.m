function [Hmo_cut_spectrum,Hmo_all,Tp]=load_FRFinst(fname,TpFlag,forecast_num_FRF)

time=ncread(fname,'time');
Hs=ncread(fname,'waveHs');
if TpFlag==2
    Tp=ncread(fname,'waveTpPeak');
    f=ncread(fname,'waveFrequency');
    theta=0;
    E_D_all=ncread(fname,'waveEnergyDensity');
    E_D_all(:,:,1)=E_D_all;
else
    Tp=ncread(fname,'waveTp');
    f=ncread(fname,'waveFrequency');
    theta=ncread(fname,'waveDirectionBins');
    E_D_all=ncread(fname,'directionalWaveEnergyDensity');
end

min_period=6;  % min allowable period
min_theta=0; % min allowable theta
max_theta=160; % max allowable theta

Nt=find(time>=forecast_num_FRF,1);
E_D=squeeze(E_D_all(:,:,Nt));

    time_reference = datenum('1970', 'yyyy');
        time_EDT = time_reference + double(time(Nt))/24/60/60-5/24;  % EDT time
        
        str1=['Nowcast Time: ' datestr(time_EDT,'yyyy-mm-dd HH:MM') ' EDT']

f_max=1/min_period;
f_max_ind=find(f>f_max,1);
t_min_ind=find(theta>min_theta,1);
t_max_ind=find(theta>max_theta,1);

f=f(1:f_max_ind);
theta=theta(t_min_ind:t_max_ind);
if TpFlag==2
    theta=0;
    t_min_ind=1;
    t_max_ind=1;
end
E_D=E_D(1:f_max_ind,t_min_ind:t_max_ind);

nf=length(f);

Hmo=0;
for i=1:nf
   for j=1:length(theta)
      if i==1
         del_f=(f(2)-f(1));
      elseif i==length(f)
         del_f=(f(length(f))-f(length(f)-1));
      else
         del_f=(f(i+1)-f(i))/2.+(f(i)-f(i-1))/2.;
      end
      
      if j==1
         if length(theta)==1
             del_theta=1;
         else
             del_theta=(theta(2)-theta(1));
         end
      elseif j==length(theta)
         del_theta=(theta(length(theta))-theta(length(theta)-1));
      else
         del_theta=(theta(j+1)-theta(j))/2.+(theta(j)-theta(j-1))/2.;
      end
      
      Hmo = Hmo + E_D(i,j)*del_f*del_theta;
   end
end
Hmo_cut_spectrum = sqrt(Hmo)*4.004;
Hmo_all=Hs(Nt);
Tp=Tp(Nt);


