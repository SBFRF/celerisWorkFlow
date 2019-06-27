function load_nc_duck_CMTB(dateIn, method)
% date in expected in epoch times
%  INPUT: 
%       dateIN is a epoch time (time in seconds since 1970-01-01)
%       method: binary, 1 is closest in time 
%                       0 is closest in history
%
% will go to server and pull bathy Closest in HISTORY 
% to the date in, then make it an appropriate bathy 
% for celeris simulations, it also sets the Gauge locations 
% 
% written by Patrick Lynett, USC modified by Spicer Bak, USACE 
%%
xMax = 910; % 
yMax = 1400; %smooth beyond this

clf
time_reference = datenum('1970', 'yyyy');
url='https://chlthredds.erdc.dren.mil/thredds/dodsC/cmtb/integratedBathyProduct/survey/survey.ncml';
    
disp([' BATHY: Loading bathy datafile: ' url ])
x=ncread(url,'xFRF', 1, inf);
y=ncread(url,'yFRF', 1, inf );
time_all=ncread(url,'time', 1, inf);  % in epoch 
ye=find(y>yMax,1);  % find portions of data to remove 
xe=find(x>xMax,1);  % find portions of data to remove 
cd ..  % necessary due to needing to be in this folder to see this script,

if method == 1
    % closest in time 
    time_diff = abs(time_all-dateIn);
    min_time_diff = min(time_diff);
    ind_min_time_diff=find(time_diff == min_time_diff,1);
elseif method == 0
    % closest in history 
    time_diff = dateIn - time_all;
    max_time_diff = min(time_diff(time_diff > 0));
    ind_min_time_diff = find(time_diff == max_time_diff,1); 
end 

time_matlab = time_reference + double(time_all(ind_min_time_diff))/24/60/60; 
cur_date=datestr(time_matlab,'yyyy-mm-ddTHH.MM.SS.Z');
disp([' BATHY: Using CMTB bathy from ' cur_date ' for simulation Date ' datestr(datetime(dateIn, 'convertfrom','posixtime')) ])
% only pull one bathy 
hAll=ncread(url,'elevation', [1, 1, ind_min_time_diff], [xe, ye, 1]);
h=squeeze(hAll)';   
disp([' BATHY: bathy Gathered, now Preprocessing'])
[ny,nx]=size(h);

% remove data now 
x=x(1:xe);
y=y(1:ye);
h=h(1:ye,1:xe);

for i=1:length(x)
    for j=1:length(y)
        if h(j,i)<-8.7
            h(j,i)=-8.7;
        end
    end
end

instruments=[900,935; 605, 938; 400, 940; 300, 940; 200, 940; 150, 940; 125, 940; 100, 940; 75, 940; 55, 940; x(1), 940];
% 8m-awac, 6m-awac, 4.5m-awac, 3.5m-adop, 2.8m-pres, 2.1m-pres, 1.9m-pres

%********************************************************************
%% plot the data
figure(1)
clf
pcolor(x,y,h)
view(0,90)
shading interp
grid off
axis equal
colorbar
xlabel('X (m)')
ylabel('y (m)')
hold on
plot(instruments(:,1),instruments(:,2),'ro')
pause(.1)

%*******************************************************************


% save for Celeris
xoffset=x(1);
yoffset=y(1);

x=x'-xoffset;
y=y'-yoffset;
B=h';

instruments(:,1)=instruments(:,1)-xoffset;
instruments(:,2)=instruments(:,2)-yoffset;

save celeris_bathy.mat x y B instruments xoffset yoffset








