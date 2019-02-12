break_loop=0;

while break_loop==0
    
    clear all;
    close all;
    w = warning ('off','all');
    clickableMap = false;
    
    break_loop=0;  % run forever
    
    %INPUTS:
    
    % default values for all simulations:
    Courant=0.1;  % Courant number, controls time step
    Mannings_n=0.005;  % Friction factor, quadratic law friction factor f/2*H*u^2
    min_depth=1e-3; % minimum depth allowable, in m, typically < 1/1000 of "offshore depth"
    
    % Forecast time
    forecast_date=['2019-01-14T13.00.00.Z']; % hindcast date in Z
    
    % screen-size
    ss=get(0,'ScreenSize');
    % resol=[1920 830];
    % resol=[2560 1600];
    resol=[ss(3) ss(4)];  % automatically detect screen resolution
    
    % wall clock time to run
    sim_time=60*60; % wall clock time
    num_frames=30; % num frames to inlcude in animation
    
    % END TOP INPUT BLOCK
    %%
    % Database locations and corresponding "ID"
    numgrids=1;
    dir_names=cell(numgrids,1);
    wave_bc=zeros(numgrids,1);
    grid_size=wave_bc;
    tidestatID=dir_names;
    
    % grid info
    ind_c=1;
    dir_names{ind_c}='Duck_FRF';
    wave_bc(ind_c)=2;
    grid_size(ind_c)=1.; % this is the numerical grid resolution in m
    tidestatID{ind_c}='8651370'; %locstr='Duck_FRF';
    
    
    % Decide on location
    % from bathy database, "choice" varaiable corresponds to "ID" in "labels"
    % variable
    
    cd_home=cd;
    cd_http=[cd '\output_http\'];
    cur_timeZ=now+datenum(2010,1,1,7,0,0)-datenum(2010,1,1,0,0,0);  % simulation / forecast time
    
    %     % update forecast file
    %     try
    %         disp(['Downloading CDIP Forecast Files'])
    %         urlwrite('http://thredds.cdip.ucsd.edu/thredds/fileServer/cdip/model/MOP_alongshore/socal_alongshore_forecast.nc','socal_alongshore_forecast.nc');
    %         urlwrite('http://thredds.cdip.ucsd.edu/thredds/fileServer/cdip/model/MOP_alongshore/norcal_alongshore_forecast.nc','norcal_alongshore_forecast.nc');
    %     end
    %
    
    disp(['Loading Input Data and Starting Celeris...'])
    
    % Load location-specifc parameters
    choice=1;
    cd(['bathy\' char(dir_names(choice))])  % change to local database directory, where the "celeris_bathy.mat" is located
    
    % default values, will change only wave input boundary
    bc(1)=1; %westBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
    bc(2)=2; %eastBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
    bc(3)=2; %southBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
    bc(4)=2; %northBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
    sponge(1)=2; % westBoundary sponge width, only used if bc(1)=2
    sponge(2)=30; % eastBoundary sponge width, only used if bc(2)=2
    sponge(3)=30; % southBoundary sponge width, only used if bc(3)=2
    sponge(4)=30; % northBoundary sponge width, only used if bc(4)=2
    
    wave_boundary_str=num2str(wave_bc(choice));
    
    for ii=1:length(wave_boundary_str)
        
        wave_boundary_c=str2num(wave_boundary_str(ii));
        
        if wave_boundary_c==1 % waves through west boundary
            bc(1)=4; %westBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
        end
        
        if wave_boundary_c==2 % waves through east boundary
            bc(2)=4; %eastBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
        end
        
        if wave_boundary_c==3 % waves through south boundary
            bc(3)=4; %southBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
        end
        
        if wave_boundary_c==4 % waves through north boundary
            bc(4)=4; %northBoundary, value=1: Solid wall; value=2: Sponge layer; value=3: Sine wave through this boundary; value=4: Irregular wave through this boundary
        end
        
    end
    
    % above theta values can also be overwritten if directions are known
    
    
    % load the pre-formatting mat file containing bathy/topo data
    load celeris_bathy.mat
    
    cd(cd_home)
    %------
    
    % load FRF THREDDS data, need coords from celeris_bathy.mat
    
    fname_year=forecast_date(1:4);
    fname_month=forecast_date(6:7);
    if length(fname_month)==1
        fname_month=['0' fname_month];
    end
    
    fname_thredds=['FRF-ocean_waves_8m-array_' fname_year fname_month '.nc'];
    fname_writeout='FRF_8m-array.nc';
    eval([' ! wget --output-document ' fname_writeout ' --no-check-certificate "https://chlthredds.erdc.dren.mil/thredds/fileServer/frf/oceanography/waves/8m-array/' fname_year '/' fname_thredds '"'])

    % get FRF data % 8m-awac, 6m-awac, 4.5m-awac, 3.5m-adop, 2.8m-pres, 2.1m-pres, 1.9m-pres
    fname_thredds8=['FRF-ocean_waves_8m-array_' fname_year fname_month '.nc'];
    fname_writeout8='FRF_8m-array.nc';
    %eval([' ! wget --output-document ' fname_writeout8 ' --no-check-certificate "https://chlthredds.erdc.dren.mil/thredds/fileServer/frf/oceanography/waves/8m-array/' fname_year '/' fname_thredds8 '"'])
    
    fname_thredds45=['FRF-ocean_waves_awac-4.5m_' fname_year fname_month '.nc'];
    fname_writeout45='FRF_awac-45m.nc';
    eval([' ! wget --output-document ' fname_writeout45 ' --no-check-certificate "https://chlthredds.erdc.dren.mil/thredds/fileServer/frf/oceanography/waves/awac-4.5m/' fname_year '/' fname_thredds45 '"'])
    
    fname_thredds35=['FRF-ocean_waves_adop-3.5m_' fname_year fname_month '.nc'];
    fname_writeout35='FRF_adop-35m.nc';
    eval([' ! wget --output-document ' fname_writeout35 ' --no-check-certificate "https://chlthredds.erdc.dren.mil/thredds/fileServer/frf/oceanography/waves/adop-3.5m/' fname_year '/' fname_thredds35 '"'])
    
    fname_thredds19=['FRF-ocean_waves_xp125m_' fname_year fname_month '.nc'];
    fname_writeout19='FRF_xp-19m.nc';
    eval([' ! wget --output-document ' fname_writeout19 ' --no-check-certificate "https://chlthredds.erdc.dren.mil/thredds/fileServer/frf/oceanography/waves/xp125m/' fname_year '/' fname_thredds19 '"'])
    
    % get webcam image at this time
    eval(' ! wget --output-document webcam_FRF.jpg "http://www.frf.usace.army.mil/oscar/nowc4.jpg"')
    webcam_image=imread('webcam_FRF.jpg');
    [wcx,wcy,wcz]=size(webcam_image);
    wc_reduce=7;
    webcam_image=webcam_image(1:wc_reduce:wcx,1:wc_reduce:wcy,:);
    [wcx,wcy,wcz]=size(webcam_image);
    pier_image=imread('pier.jpg');
    [px,py,pz]=size(pier_image);
    
    
    load_FRFwave(fname_writeout)
    load FRFwave_forecast.mat
    
    disp(['Water depth at FRF nowcast offshore boundary (m): ' num2str(nominaldepth)])
    
    time_reference = datenum('1970', 'yyyy');
    
    % find FRFwave data indices that we will generate forecasts for
    forecast_count=1;  % only forecast for most recent data, can loop across this variable if more times desired
    cur_dateZ=datestr(cur_timeZ,'yyyy-mm-ddTHH.MM.SS.Z');
    forecast_num_matlab=datenum(forecast_date,'yyyy-mm-ddTHH.MM.SS.Z');
    forecast_num_FRF=(forecast_num_matlab-time_reference)*24*60*60;
    Nt=find(waveTime>=forecast_num_FRF,1);
    if isempty(Nt)==1  % if forecast_time is greater than last waveTime, model for last waveTime
        Nt=length(waveTime);
    end
    
    time_matlab = time_reference + double(waveTime(Nt))/24/60/60;  % for some reason time stamps off by 4 min 15 sec
    cur_date=datestr(time_matlab,'yyyy-mm-ddTHH.MM.SS.Z');
    FRF_Wt_ind(forecast_count)=Nt;
    
    %load local tide datay
    statID=char(tidestatID(choice));
    
    % get last few days of tide predictions and data:
%     days_past=addtodate(cur_timeZ,-3,'day');
%     days_future=addtodate(cur_timeZ,10,'day');
%     days_past_str=datestr(days_past,'yyyymmdd HH:MM');
%     days_now_str=datestr(cur_timeZ,'yyyymmdd HH:MM');
%     days_future_str=datestr(days_future,'yyyymmdd HH:MM');
%     eval([' ! wget --output-document predicted_NAVD88.csv "http://tidesandcurrents.noaa.gov/api/datagetter?begin_date=' days_past_str '&end_date=' days_future_str '&station=' statID '&product=predictions&datum=NAVD&units=metric&time_zone=gmt&application=USC&format=csv"'])
%     eval([' ! wget --output-document predicted.csv "http://tidesandcurrents.noaa.gov/api/datagetter?begin_date=' days_past_str '&end_date=' days_future_str '&station=' statID '&product=predictions&datum=mllw&units=metric&time_zone=gmt&application=USC&format=csv"'])
%     eval([' ! wget --output-document measured_NAVD88.csv "http://tidesandcurrents.noaa.gov/api/datagetter?begin_date=' days_past_str '&end_date=' days_now_str '&station=' statID '&product=one_minute_water_level&datum=NAVD&units=metric&time_zone=gmt&application=USC&format=csv"'])
%     eval([' ! wget --output-document measured.csv "http://tidesandcurrents.noaa.gov/api/datagetter?begin_date=' days_past_str '&end_date=' days_now_str '&station=' statID '&product=one_minute_water_level&datum=mllw&units=metric&time_zone=gmt&application=USC&format=csv"'])

    % FRF tide data
    fname_thredds_tide=['FRF-ocean_waterlevel_eopNoaaTide_' fname_year fname_month '.nc'];
    fname_writeout_tide='FRF_tides.nc';
    eval([' ! wget --output-document ' fname_writeout_tide ' --no-check-certificate "https://chlthredds.erdc.dren.mil/thredds/fileServer/frf/oceanography/waterlevel/eopNoaaTide/' fname_year '/' fname_thredds_tide '"'])

    tideTime=ncread(fname_writeout_tide,'time');
    waterlevel=ncread(fname_writeout_tide,'waterLevel');
    pred_waterlevel=ncread(fname_writeout_tide,'predictedWaterLevel');
    
%     % load predictions
%     TC_pred=importdata('predicted.csv',',',1);
%     tide_pred=TC_pred.data(:,1);
%     m=length(tide_pred);
%     tide_pred=tide_pred(2:m);
%     date_pred=datenum(TC_pred.textdata(2:m,1),'yyyy-mm-dd HH:MM');
%     
    % load predictions NAVD88
    tide_pred=pred_waterlevel;
%     m=length(tide_pred);
%     tide_pred=tide_pred(2:m);
    
    % load observations
 %   TC_meas=importdata('measured.csv',',',1);
 %   tide_meas=TC_meas.data(:,1);
 %   m=length(tide_meas);
  %  tide_meas=tide_meas(2:m);
   % date_meas=datenum(TC_meas.textdata(2:m,1),'yyyy-mm-dd HH:MM');
    
    % load observations NAVD88
    tide_meas=waterlevel;
%     m=length(tide_meas);
%     tide_meas=tide_meas(2:m);
    
    %%%%%
    % add value here for storm surge, tide - any sort of constant water level
    % change for entire domain
    water_level_change=0;  %(m), all grids at MHW datum, shift down 0.4' to MSL for surge inclusion, than add surge, dummy value here for now
    load_bathy  % run bathy load script
    bathyimage=imread('bathytopo.jpg');
    [bx,by,bz]=size(bathyimage);
    
    
    %clean and create output directory
    grid_name=char(dir_names(choice));
    frames_dir=[cd '\output\' grid_name];
    try
        eval(['rmdir ' frames_dir ' s;'])
    end
    pause(2)  % allow dropbox to figure it out
    mkdir(frames_dir)
    
    % start forecast loop
    for ifcst=1:forecast_count
        Nt=FRF_Wt_ind(ifcst);
        time_matlab = time_reference + double(waveTime(Nt))/24/60/60;  % for some reason time stamps off by 4 min 15 sec
        cur_date=datestr(time_matlab,'yyyy-mm-ddTHH.MM.SS.Z');
        
        % create spectrum file
        hf3=figure(3);
        H=Hs(Nt);
        T=Tp(Nt);
        theta=Dp(Nt);
        
        spectrum_FRF_2D_interp(waveFrequency,waveMeanDirection,squeeze(waveEnergyDensity(:,:,Nt))',H,T)
        
        set(hf3,'PaperPosition',[0 0 3 3]*resol(1)/2560);
        print -djpeg100 spectrum2D.jpg
        specimage=imread('spectrum2D.jpg');
        [wx,wy,wz]=size(specimage);
        
        % create forecast plot
        hf5=figure(5);
        clf
        is=1;
        ie=length(waveTime);
        wT_trunc=time_reference + double(waveTime(is:ie))/24/60/60;
        meanH=mean(Hs(is:ie));
        
        
        mint_plot=forecast_num_matlab-3.5;
        maxt_plot=forecast_num_matlab+12/24 ;
        
        subplot(4,1,1)
        hold on
        plot([cur_timeZ cur_timeZ],[-10 50],'--g','LineWidth',1)
        plot(time_reference + double(waveTime(Nt))/24/60/60,Hs(Nt),'r.','MarkerSize',15)
        plot(wT_trunc,Hs(is:ie))
        axis([mint_plot maxt_plot 0 1.1*max(Hs(is:ie))])
        datetick('x','mm-dd-HH','keepticks')
        grid on
        ylabel('Significant Wave Height (m)','FontSize',5,'fontweight','bold')
        xlabel('Date/Time (UTC)','fontsize',5,'fontweight','bold')
        
        time_EDT = time_reference + double(waveTime(Nt))/24/60/60-5/24;  % EDT time
        
        str1=['Nowcast Time: ' datestr(time_EDT,'yyyy-mm-dd HH:MM') ' EDT'];
        str2=['Plot Generated on: ' datestr(now,'yyyy-mm-dd HH:MM') ' PDT'];
        
        title({str1; str2},'FontSize',5,'fontweight','bold')
        set(gca,'fontsize',5)
        
        subplot(4,1,2)
        hold on
        plot([cur_timeZ cur_timeZ],[-10 50],'--g','LineWidth',1)
        plot(time_reference + double(waveTime(Nt))/24/60/60,Tp(Nt),'r.','MarkerSize',15)
        plot(wT_trunc,Tp(is:ie))
        axis([mint_plot maxt_plot 0 1.1*max(Tp(is:ie))])
        datetick('x','mm-dd-HH','keepticks')
        grid on
        ylabel('Peak Wave Period (sec)','FontSize',5,'fontweight','bold')
        xlabel('Date/Time (UTC)','fontsize',5,'fontweight','bold')
        legend('Current Time','Animation Time','Location','SouthWest')
        set(gca,'fontsize',5)
        
        subplot(4,1,3)
        hold on
        plot([cur_timeZ cur_timeZ],[-10 500],'--g','LineWidth',1)
        plot(time_reference + double(waveTime(Nt))/24/60/60,Dp(Nt),'r.','MarkerSize',15)
        plot(wT_trunc,Dp(is:ie))
        axis([mint_plot maxt_plot 0.99*min(Dp(is:ie)) 1.01*max(Dp(is:ie))])
        datetick('x','mm-dd-HH','keepticks')
        grid on
        ylabel('Wave Direction (deg CW from N)','FontSize',5,'fontweight','bold')
        xlabel('Date/Time (UTC)','fontsize',5,'fontweight','bold')
        set(gca,'fontsize',5)
        
        subplot(4,1,4)
        hold on
        plot([cur_timeZ cur_timeZ],[-10 500],'--g','LineWidth',1)
        
        water_level_change=interp1(tideTime,tide_meas,waveTime(Nt));
        
        plot(time_reference + double(waveTime(Nt))/24/60/60,water_level_change,'r.','MarkerSize',15)
        
        plot(time_reference + double(tideTime)/24/60/60,tide_pred)
        plot(time_reference + double(tideTime)/24/60/60,tide_meas,'k')
        axis([mint_plot maxt_plot 1.1*min(tide_pred) 1.1*max(tide_pred)])
        datetick('x','mm-dd-HH','keepticks')
        grid on
        ylabel('Tide (m - NAVD88)','FontSize',5,'fontweight','bold')
        xlabel('Date/Time (UTC)','fontsize',5,'fontweight','bold')
        %        legend('Current Time','Observed Tides','Predicted Tides','Location','SouthEast')
        set(gca,'fontsize',5)
        
        set(hf5,'PaperPosition',[0 0 3.2 7.1]*resol(1)/2560*0.9194);
        print -djpeg100 FRF_forecast.jpg
        forecastimage=imread('FRF_forecast.jpg');
        [fcx,fcy,fcz]=size(forecastimage);
        % trimoff top of image
        forecastimage=forecastimage(15:fcx-20,:,:);
        [fcx,fcy,fcz]=size(forecastimage);
        
        % write  for web page
        set(hf5,'PaperPosition',[0 0 3.2 7.1]*resol(1)/2560*0.9194*1);
        print -djpeg100 FRF_forecast_small.jpg
        forecastimage_small=imread('FRF_forecast_small.jpg');
        [fcx_s,fcy_s,fcz_s]=size(forecastimage_small);
        forecastimage_small=forecastimage_small(20:fcx_s-40,:,:);
        fnameforecast_c=[grid_name '_Celeris_' cur_date '_forecast.jpg'];
        cd(frames_dir)
        imwrite(forecastimage_small,fnameforecast_c,'jpg','Quality',75);
        cd(cd_home)
        
        % regenerate bathy to account for tide level
        disp(['Shifting bathy to account for tide level of (m-NAVD88): ' num2str(water_level_change)])
        cd(['bathy\' char(dir_names(choice))])  % change to local database directory, where the "celeris_bathy.mat" is located
        load celeris_bathy.mat
        cd(cd_home)
        
        load_bathy  % run bathy load script
        bathyimage=imread('bathytopo.jpg');
        [bx,by,bz]=size(bathyimage);
        
        % find instruement node indices for output
        inst_ind=instruments*0;
        for i=1:length(instruments(:,1))
            inst_ind(i,1)=find(x>=instruments(i,1),1);
            inst_ind(i,2)=find(y>=instruments(i,2),1);
        end
        
        % strip range for output
        range_x=[min(inst_ind(:,1)) max(inst_ind(:,1))];
        range_y=median(inst_ind(:,2));
        
        % Determine time step
        max_depth=-min(min(B));
        ds=min(mean(diff(x)),mean(diff(y)));
        timestep=Courant*ds/sqrt(9.81*max_depth);
        
        % write the simulation control file
        file_name_cml='matlab_launch.cml';
        
        % create colorbar for image overlay
        % cap color red limit at 2.5 m
        H_plot_cbar=min(2.5/2,meanH);
        colorbar_image=create_colorbar_jpg([-2.0*H_plot_cbar/2 4.0*H_plot_cbar/2]);
        colorbar_image=imrotate(colorbar_image,-90);
        [cbx,cby,cbz]=size(colorbar_image);
        
        % set breaking parameters based on wave height.  These are VIZ
        % not pyhsical model parameters - they dont change simulation
        % results
        brk_limH=[0.65 1.5];
        brk_limthres=[0.2 0.45];
        plung_atten=0.02;
        
        if H>brk_limH(2) % energetic plunging
            brk1=brk_limthres(2);  % on/off threshold
        elseif H<brk_limH(2) % low energy spilling
            brk1=brk_limthres(1);
        else  % seomthing inbetween
            brk1=interp1(brk_limH,brk_limthres,H);
        end
        brk2=brk_limthres(2)-brk1+plung_atten;  % attenuation
        
        
        %------------
        % Type of input wave !not used for coastal spectrum-driven
        % waves!. just putting dummy values
        wave_type=2; %=1 for solitary wave, =2 for sine wave through boundary
        xco=0; % defines the initial x-location (m) of the crest of the solitary wave
        yco=0; % defines the initial y-location (m) of the crest of the solitary wave
        
        create_cml(file_name_cml,timestep,Mannings_n,m_width,m_length,nx,ny,wave_type,H_plot_cbar,T,theta,xco,yco,bc,sponge,min_depth,brk1,brk2,inst_ind,range_x,range_y)
        
        delete array.txt
        delete time_axis.txt
        
        % run Celeris
        %           ! Celeris.exe  &
        
        % for some reason the straight call doesnt always maximize.
        % Force with robot
        
        import java.awt.*;
        import java.awt.event.*;
        rob=Robot;
        
        h = actxserver('WScript.Shell');
        h.Run('Celeris.exe'); %Invokes
        
        pause_int=sim_time/2/num_frames;
        pause(sim_time/2+30)  % wait for sim warm up
        
        h.AppActivate('Celeris.exe'); %Brings to focus
        rob.keyPress(KeyEvent.VK_WINDOWS)
        rob.keyPress(KeyEvent.VK_UP)
        rob.keyRelease(KeyEvent.VK_WINDOWS)
        rob.keyRelease(KeyEvent.VK_UP)
        
        pause(10)
        
        logoimage=imread('logo.jpg');
        [lx,ly,lz]=size(logoimage);
        
        
        cd_home=cd;
        for i=1:num_frames
            imageData = screencapture(0, [0,0,resol(1)-500,resol(2)]); % select a small desktop region, crop out user menu
            [nxiD,nyiD,nziD]=size(imageData);
            
            if i==1
                imageData_anim=zeros(nxiD,nyiD,nziD,num_frames,'uint8');
            end
            imageData(441:440+px,371:370+py,:)=pier_image;
            imageData=imrotate(imageData, 18.2,'crop'); 
            
            %           imageData(nxiD-bx+1:nxiD,1:by,:)=bathyimage;
            imageData(nxiD-wx+1:nxiD,nyiD-wy+1:nyiD,:)=specimage;
            imageData(1:cbx,nyiD-cby+1-300:nyiD-300,:)=colorbar_image;
            imageData(cbx+21:wcx+cbx+20,nyiD-wcy+1:nyiD,:)=webcam_image;
            imageData(1:fcx,1:fcy,:)=forecastimage;
            imageData(nxiD-lx+1:nxiD,1:ly,:)=logoimage;
            imageData_anim(:,:,:,i)=imageData;  % straight screen capture stack
            
            pause(pause_int)
        end
        
        ! taskkill /IM Celeris.exe > screen.txt
        
        % write final surface to image
        cd(frames_dir)
        fname_c=[grid_name '_Celeris_' cur_date '.jpg'];
        imwrite(imageData,fname_c); % save the captured image to file
        
        % create animated gif
        anim_name=[grid_name '_Celeris_' cur_date '.gif'];
        for i=1:num_frames
            im = imageData_anim(:,:,:,i);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if i == 1
                imwrite(imind,cm,anim_name,'gif','DelayTime',0,'Loopcount',inf);
            else
                imwrite(imind,cm,anim_name,'gif','WriteMode','append','DelayTime',0);
            end
        end
        
        % create animated thumbnail gif
        anim_name=[grid_name '_Celeris_' cur_date '_thumb.gif'];
        reduce_resfac=4;
        i_crop=100;
        j_crop=340;
        for i=max(1,num_frames-80):4:num_frames
            im = imageData_anim(i_crop:reduce_resfac:nxiD-i_crop,j_crop:reduce_resfac:nyiD-j_crop,:,i);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if i == max(1,num_frames-80)
                imwrite(imind,cm,anim_name,'gif','DelayTime',0,'Loopcount',inf);
            else
                imwrite(imind,cm,anim_name,'gif','WriteMode','append','DelayTime',0.3);
            end
        end
        
        % write animation to file
        anim_name=[grid_name '_Celeris_' cur_date];
        outputVideo = VideoWriter(fullfile(cd,anim_name));
        outputVideo.FrameRate = 10;
        outputVideo.Quality = 30;
        open(outputVideo);
        writeVideo(outputVideo,imageData_anim);
        close(outputVideo);
        cd(cd_home)
        clear imageData_anim;
        clear imageData_scap;
    end
    
    cd(frames_dir)
    fid = fopen('model_index.html','w');
    fprintf(fid, ['<HTML><HEAD><TITLE>USC Coastal Wave Nowcast - ' grid_name '</TITLE>\n']);
    %fprintf(fid, ['<link rel="stylesheet" href="exp.css" type="text/css">\n']);
    fprintf(fid, ['</head>\n']);
    %fprintf(fid, ['<BODY background="../../bg.gif">\n?);
    %fprintf(fid, ['<A NAME="top" HREF="../index.html">CCL - POLB Wave Prediction Home Page</A>\n']);
    %fprintf(fid, ['<H1>CDIP 215 Data Inventory - Last 30 Days</H1>\n']);
    %fprintf(fid, ['<HR SIZE=5>\n']);
    
    fprintf(fid, ['<BR> ---------- ' ]);
    fprintf(fid, ['<BR> ' grid_name ' Wave Nowcast ']);
    fprintf(fid, ['<BR> Last updated at ' datestr(now,'yyyy-mm-dd HH:MM') ' PDT ' ]);
    fprintf(fid, ['<BR> Note: If image links appear broken (images not showing) and the last update time is recent, files are currently being replaced. Try again in a few minutes. ' ]);
    fprintf(fid, ['<BR> ---------- ' ]);
    
    for ifcst=1:forecast_count
        Nt=FRF_Wt_ind(ifcst);
        time_matlab = time_reference + double(waveTime(Nt))/24/60/60;  % for some reason time stamps off by 4 min 15 sec
        cur_date=datestr(time_matlab,'yyyy-mm-ddTHH.MM.SS.Z');
        
        time_EDT = time_matlab-5/24;  % for some reason time stamps off by 4 min 15 sec
        cur_date_EDT=datestr(time_EDT,'yyyy-mm-ddTHH.MM.SS');
        
        fname_jpg=[grid_name '_Celeris_' cur_date '.jpg'];
        fname_gif=[grid_name '_Celeris_' cur_date '.gif'];
        fname_gifthumb=[grid_name '_Celeris_' cur_date '_thumb.gif'];
        anim_name=[grid_name '_Celeris_' cur_date '.avi'];
        fnameforecast_c=[grid_name '_Celeris_' cur_date '_forecast.jpg'];
        image_width=700;
        
        fprintf(fid, ['<BR> <a href="./' fname_gif '"> <img src="./' fnameforecast_c '" height="' num2str(nxiD/nyiD*image_width)  '" width="' num2str(nxiD/nyiD*image_width*fcy/fcx)  '" </a>> <a href="./' fname_gif '"> <img src="./' fname_gifthumb '" height="' num2str(nxiD/nyiD*image_width)  '" width="' num2str(image_width)  '"> </a>']);
        fprintf(fid, ['<BR> Click on image above to view wave field animation (animated gif), or download avi file below']);
        fprintf(fid, ['<BR> <a href="./' fname_jpg '">' fname_jpg '</a> -- Wave Surface Snapshot at forecast time ' cur_date_EDT ' EDT']);
        fprintf(fid, ['<BR> <a href="./' anim_name '">' anim_name '</a> -- Wave Surface Animation at forecast time ' cur_date_EDT ' EDT']);
        fprintf(fid, ['<BR> ---------- ' ]);
    end
    fclose(fid);
    cd(cd_home)
    
    % copy files to http
    http_dir=[cd_http grid_name];
    eval(['! del /s/q ' http_dir '\* > screen_copy.txt'])
    pause(2);
    eval(['copyfile ' frames_dir '\model_index.html ' http_dir])
    eval(['copyfile ' frames_dir '\*.* ' http_dir])
    
    break_loop=1;
    
    load_array
end






