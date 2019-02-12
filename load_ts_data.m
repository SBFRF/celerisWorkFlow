function load_ts_data


start_analysis_time=1800; % time series analysis start time, in seconds


load ts_locations.dat
[m_ts,n_ts]=size(ts_locations);     
cd tmsr
for loc=1:m_ts
   ind=['000' num2str(loc)];
   file_ind=ind(length(ind)-3:length(ind));
   
   filename=['tmsr' file_ind '.dat'];
   count=1;
   while exist(filename)==0
      ind=['000' num2str(loc-count)];
      file_ind=ind(length(ind)-3:length(ind));
      
      filename=['tmsr' file_ind '.dat']; 
      count=count+1;
   end
   
   data=0;
   eval(['load tmsr' file_ind '.dat'])
   eval(['data=tmsr' file_ind ';'])
   eval(['clear tmsr' file_ind ';'])
   [m,n]=size(data);
   
   
   xloc(loc)=data(1,1);
   yloc(loc)=data(1,2);
   hloc(loc)=data(1,3);
   
   t_s=data(2:m-1,1);
   z_s=data(2:m-1,2);
   
   dt_s=t_s(5)-t_s(4);
   t_s=t_s-start_analysis_time;
   time=[0:dt_s:t_s(length(t_s))];
   eta=interp1(t_s,z_s,time);
   [Hs(loc),zmean(loc),Tp(loc),Tm(loc),f,Se,f_ave,S_ave]=spectral(time,eta);   
   
   [loc,m_ts,Hs(loc)]
   pause(0.0001)
end
cd ..
save output.mat



