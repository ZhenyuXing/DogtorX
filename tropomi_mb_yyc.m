%% 
% The following section is to read the data from raw files (path). Valid observations 
% allocated in the predefined grid cells are averaged for each individual file 
% (path). The mean_ch4 is the temporal averaged XCH4 based on multiple files (paths) 
% during the predefined period. The coverage_count means the total number of valid 
% observations in each grid cell when do the temporal averaging.
% 
% Data reading is designed for year-by-year processing! Repeat for 2021 and 
% 2022 to obtain 3-year observations. 

% For user input: year, start month, start date, end month, end date, etc., file path.
clear; clc;
prompt={'Enter the FOUR digit year-yyyy:',...
    'Enter the TWO digit start month-mm:','Enter the TWO digit start date:',...
    'Enter the TWO digit end month:','Enter the TWO digit end date:',...
    'Enter the minimum latitude of interest:','Enter the maximum latitude of interest:',...
    'Enter the minimum longitude of interest:','Enter the maximum longitude of interest:',...
    'Enter the qa threshold value (0 to 1 indicates worst to best):'};
name='For user input';
numlines=1;
defaultanswer={'2020','01','01','12','31','50.0','52.0','-115.0','-113.0','0.5'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
u_ipt = inputdlg(prompt,name,numlines,defaultanswer,options);
yyyy = str2double(u_ipt{1});startmm = str2double(u_ipt{2});startdd = str2double(u_ipt{3});
endmm = str2double(u_ipt{4});enddd = str2double(u_ipt{5});minlat = str2double(u_ipt{6});
maxlat = str2double(u_ipt{7});minlon = str2double(u_ipt{8});maxlon = str2double(u_ipt{9});
qa_pass = str2double(u_ipt{10});
folder_name = uigetdir('','Could you please tell me whare are your .nc files?');
flist = getfilelist('CH4',folder_name,yyyy,startmm,startdd,endmm,enddd); %Get .nc files list associated with the user input time-period
[grid_lon,grid_lat,ch4,mean_ch4,coverage_count,flist_ex] = load_ch4(flist, minlat, maxlat, minlon, maxlon, qa_pass);
save('D:\Urban_methane\data\processed_data\mat\yyc_2020.mat'); % User-defined path.
%% This section is to integrate the 3 year observations.
cc_yyc_2020(isnan(cc_yyc_2020)) = 0;
cc_yyc_2021(isnan(cc_yyc_2021)) = 0;
cc_yyc_2022(isnan(cc_yyc_2022)) = 0;
cc_yyc = cc_yyc_2020 + cc_yyc_2021 + cc_yyc_2022;
ch4_yyc = ch4_yyc_2020;
ch4_yyc(:,:,size(ch4_yyc,3)+1:size(ch4_yyc,3)+size(ch4_yyc_2021,3)) = ch4_yyc_2021;
ch4_yyc(:,:,size(ch4_yyc,3)+1:size(ch4_yyc,3)+size(ch4_yyc_2022,3)) = ch4_yyc_2022;
%% This section is to determine the boundaries of region of interest.
% shpefiles are provided.
countries = shaperead('D:\Emission_estimate\TROPOMI_CH4\processed_data\GIS_data\World_Countries\World_Countries__Generalized_.shp');
Canada = shaperead('D:\Emission_estimate\TROPOMI_CH4\processed_data\GIS_data\lpr_000a21a_e\lpr_000a21a_e.shp');
info = shapeinfo('D:\Emission_estimate\TROPOMI_CH4\processed_data\GIS_data\lpr_000a21a_e\lpr_000a21a_e.shp');
proj = info.CoordinateReferenceSystem;
for i = 1:13
    lat = Canada(i).X;
    lon = Canada(i).Y;
    lat(isnan(lat)) = [];
    lon(isnan(lon)) = [];
    [lat,lon] = projinv(proj,lat,lon);
    Canada(i).X = lat;
    Canada(i).Y = lon;
end
clear lat lon i
Alberta = shaperead('D:\Emission_estimate\TROPOMI_CH4\processed_data\GIS_data\lcsd000b16a_e\lcsd000b16a_e.shp');
info = shapeinfo('D:\Emission_estimate\TROPOMI_CH4\processed_data\GIS_data\lcsd000b16a_e\lcsd000b16a_e.shp');
proj = info.CoordinateReferenceSystem;
for i = 1:length(Alberta)
    if strcmp(Alberta(i).CSDNAME,'Calgary') == 1
        ax = [Alberta(i).X];ay = [Alberta(i).Y]; % Calgary lon (X) and lat (y)
        [c_lat,c_lon] = projinv(proj,ax,ay);
    else
        continue
    end
end
roi_wc = struct;
roi_wc(1).Latitude = c_lat;
roi_wc(1).Longitude = c_lon;
roi_wc(1).Region = 'Calgary';
roi_wc(2).Latitude = [50.5 50.5 51.5 51.5 50.5];
roi_wc(2).Longitude = [-114.5 -113.5 -113.5 -114.5 -114.5];
roi_wc(2).Region = 'Calgary_b';
%% The following section is to read the wind speed from the ground-based stationary observations.
% User should firstly download data from the provided data sources!

ws_yyc_ams = readtable('D:\Urban_methane\data\ams\continuous_monitoring\For_WindTrax\Multiple-Stations-wind-speed.xlsx',...
    'Range','C3:E11018','ReadVariableNames',false);
ws_yyc_ams = table2array(ws_yyc_ams);
% Only select the wind speed recorded during the time window that TROPOMI passes over YYC.
ws_mean_hourly = zeros(24,3);
ws_mean_hourly(ws_mean_hourly==0) = nan;
ws_std_hourly = zeros(24,3);
ws_std_hourly(ws_std_hourly==0) = nan;
for i = 1:24
    ws_mean_hourly(i,:) = mean(ws_yyc_ams(i:24:10993+i-1,:),1,"omitnan");
    ws_std_hourly(i,:) = std(ws_yyc_ams(i:24:10993+i-1,:),1,"omitnan");
end
ws_timewindow = [];
for i = 1:459
    ws_timewindow(end+1:end+3,:) = ws_yyc_ams(13+(i-1)*24:15+(i-1)*24,:);
end
% Valid XCH4 observations from TROPOMI are intensively distributed in warmer months. 
ws_timewindow = mean(ws_timewindow,2,"omitmissing"); % Data was recorded as km/hr
ws_mean_vc = mean(ws_timewindow,1,"omitmissing");
ws_std_vc = std(ws_timewindow,1,"omitmissing");
ws_mean_vc = ws_mean_vc*1000/60/60; % Convert km/hr to m/s;
ws_std_vc = ws_std_vc*1000/60/60; % Convert km/hr to m/s;
%% The following section is to read the surface pressure data from the MERRA-2 reanalysis data.
% Data source: https://disc.gsfc.nasa.gov/datasets/M2IMNPASM_5.12.4/summary?keywords=wind

[ wlon_20,wlat_20,uwind_mw_20,vwind_mw_20,surface_pressure_mw_20,edge_heights_mw_20 ] = ...
    merra2wind( num2str(2020),'all' );
[ wlon_21,wlat_21,uwind_mw_21,vwind_mw_21,surface_pressure_mw_21,edge_heights_mw_21 ] = ...
    merra2wind( num2str(2021),'all' );
[ wlon_22,wlat_22,uwind_mw_22,vwind_mw_22,surface_pressure_mw_22,edge_heights_mw_22 ] = ...
    merra2wind( num2str(2022),'all' );
surface_pressure_mw = ( surface_pressure_mw_22 + surface_pressure_mw_21 + surface_pressure_mw_20)/3;
wlon = wlon_21; wlat = wlat_21;
%% Applying mass balance method to calculate the emissions.
ch4_yyc_mean = mean(ch4_yyc,3,"omitmissing"); % Annual mean XCH4
ch4_yyc_std = std(ch4_yyc,0,3,"omitnan"); % Annual stdev XCH4
ch4_yyc_mean = imgaussfilt(ch4_yyc_mean,1,'FilterSize',3); % Conduct Gaussian filtering 
cc_yyc(cc_yyc==0) = nan; % replace 0 with nan in the total coverage count matrix
n_min = 10;% Threshold of valid obervation numbers for quantification, e.g. 10, 20，30.
robust_nobs = cc_yyc >= n_min;
cc_yyc = cc_yyc .* robust_nobs;
cc_yyc(cc_yyc == 0) = nan;
ch4_yyc_mean = ch4_yyc_mean .* robust_nobs;
ch4_yyc_mean(ch4_yyc_mean == 0) = nan;
ch4_yyc_std = ch4_yyc_std .* robust_nobs;
ch4_yyc_std(ch4_yyc_std == 0) = nan;
% Define the values of other parameters.
mch4 = ch4_yyc_mean;surface_pressure = surface_pressure_mw;
cell_area = 0.05 .* 111.32 .* 0.05 .* 111.32 .* cosd(grid_lat); % the pixel specific area in km^2.
M = 5.345; % M is a constant conversion factor (5.345 × 1e−9 MtCH4/km2/ppb or 5.345 kg/km2/ppb)
C = 2;
emi_mb = struct;
for i = 1
    lat = roi_wc(i).Latitude; lon = roi_wc(i).Longitude;
    lat_b = roi_wc(i+1).Latitude; lon_b = roi_wc(i+1).Longitude;
    [in,on] = inpolygon(grid_lon,grid_lat,lon,lat);
    in = in + on;
    in(in>1)=1;
    [in_b,on_b] = inpolygon(grid_lon,grid_lat,lon_b,lat_b);
    in_b = in_b + on_b;
    in_b(in_b>1)=1;
    in_sur = in;
    in_sur(in_sur==0)=2;
    in_sur(in_sur==1)=0;
    in_sur(in_sur==2)=1;
    in_sur = in_sur .* in_b;
    [in_mw,on_mw] = inpolygon(wlon,wlat,lon_b,lat_b);
    in_mw = in_mw + on_mw;
    in_mw(in_mw>1)=1;
    ch4_b = in_b .* mch4; ch4_b(ch4_b == 0) = nan;% data out of the roi filled with nan
    ch4_roi = in .* mch4; ch4_roi(ch4_roi == 0) = nan;% data out of the roi filled with nan
    area_roi = in .* cell_area;area_roi(area_roi == 0) = nan;% data out of the roi filled with nan
    mean_ws = ws_mean_vc; % the observed mean wind speed over roi, m/s.
    std_ws = ws_std_vc; % the stdev of theobserved wind speed over roi, m/s.
    surfp_roi = in_mw .* surface_pressure; surfp_roi(surfp_roi == 0) = nan;% data out of the roi filled with nan
    mean_b = mean(ch4_b(:),'omitnan');
    median_b = median(ch4_b(:),'omitnan');
    std_b = std(ch4_b(:),'omitnan');
    c = (mean_b-median_b)/std_b;
    if c>0.3
        bgd_b = median_b;
    else
        bgd_b = (2.5 * median_b) - (1.5 * mean_b);
    end
    mean_roi = mean(ch4_roi(:),'omitnan');
    std_roi = std(ch4_roi(:),'omitnan');
    delta_roi = ch4_roi - bgd_b; % the enhancement of XCH4 over roi
    delta_roi(delta_roi < 1*std_roi) = nan; % Exclude the enhancement smaller than the standard deviation.
    N = sum(~isnan(delta_roi(:))); % Obtain the number of source cells.
    d_roi(:,:,i) = delta_roi;
    std_delta_roi = std(delta_roi(:),'omitnan');
    mean_delta_roi = mean(delta_roi(:),'omitnan'); % the mean enhancement of XCH4 over roi
    M_exp = mean(surfp_roi(:),'omitnan') ./ (100*1013.0);
    L = sqrt(N*(0.05*111.32 * 0.05 * cos(mean(lat(:),'omitnan').*pi/180) * 111.32)); % the effective size of the emission region
    V = mean_ws * 60 * 60 * 24 / 1000; % convert the m/s to km/d
    V_unc = std_ws * 60 * 60 * 24 / 1000; % convert the m/s to km/d
    % The equation to calculate the emissions
    est_emi_mb = mean_delta_roi .* M .* M_exp .* L .* V .* C ./ 1000; % the unit of emission is t/d
    % The following equations to calculate the uncertainties
    est_unc_delta = std_delta_roi .* M .* M_exp .* L .* V .* C ./ 1000; % the unit is t/d
    est_unc_wind = mean_delta_roi .* M .* M_exp .* L .* V_unc .* C ./ 1000; % the unit is t/d
    est_unc_total = sqrt(est_unc_delta.^2+est_unc_wind.^2); % the unit is t/d
    % Save the parameters and emissions
    emi_mb(i).mean_xch4_b = mean_b;
    emi_mb(i).stdev_xch4_b = std_b;
    emi_mb(i).median_xch4_b = median_b;
    emi_mb(i).background_xch4 = bgd_b;
    emi_mb(i).mean_xch4_roi = mean_roi;
    emi_mb(i).stdev_xch4_roi = std_roi;
    emi_mb(i).enhancement_xch4 = mean_delta_roi;
    emi_mb(i).number_of_source_cells = N;
    emi_mb(i).area_roi = sum(area_roi(:),"all","omitnan");
    emi_mb(i).dimensionless_factor = M_exp;
    emi_mb(i).wind_speed = mean_ws;
    emi_mb(i).wind_std = std_ws;
    emi_mb(i).estimated_emission = est_emi_mb;
    emi_mb(i).estimated_uncertainty_delta = est_unc_delta;
    emi_mb(i).estimated_uncertainty_wind = est_unc_wind;
    emi_mb(i).estimated_uncertainty_total = est_unc_total;
    emi_mb(i).region = roi_wc(i).Region;
    emi_mb(i).estimated_flux = est_emi_mb*365/sum(area_roi(:),"all","omitnan"); % emission flux, t/yr/km2
    emi_mb(i).flux_uncertainty_delta = est_unc_delta*365/sum(area_roi(:),"all","omitnan");
    emi_mb(i).flux_uncertainty_wind = est_unc_wind*365/sum(area_roi(:),"all","omitnan");
    emi_mb(i).flux_uncertainty_total = est_unc_total*365/sum(area_roi(:),"all","omitnan");
end