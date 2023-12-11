function [alt_avg, alt_std] = load_alt(flist, minlat, maxlat, minlon, maxlon)
% LOAD AIR: This function is to derive column of dry air data from TROPOMI file which has been labeled with suspect plumes.
%   Input argugments should include the .nc file name and [minlat maxlat][minlon maxlon].
%   Output argument is the column of dry air in the region of interest on the corresponding day.
[grid_lon, grid_lat] = meshgrid(minlon:0.05:maxlon,minlat:0.05:maxlat);
[a,b] = size(grid_lat);
alt_avg = zeros(a,b);
alt_avg(alt_avg == 0) = nan; % replace filled 0 values with nan.
alt_std = zeros(a,b);
alt_std(alt_std == 0) = nan; % replace filled 0 values with nan.
for k = 1:length(flist)
    altitude = ncread(fname,'/PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude');
    latitude = ncread(fname,'/PRODUCT/latitude');
    longitude = ncread(fname,'/PRODUCT/longitude');
    w_l = waitbar(0,'Extracting surface altitude, please wait...');
    for j = 1:b
        for i = 1:a
            poly_lat = [minlat+0.05*(i-1) minlat+0.05*(i-1) minlat+0.05*i minlat+0.05*i minlat+0.05*(i-1)];
            poly_lon = [minlon+0.05*(j-1) minlon+0.05*j minlon+0.05*j minlon+0.05*(j-1) minlon+0.05*(j-1)];
            in = inpolygon(longitude, latitude, poly_lon, poly_lat);
            alt_raw = altitude(in);
            alt_avg(i,j) = mean(alt_raw(:),'omitnan');
            alt_std(i,j) = std(alt_raw(:),'omitnan');
        end
    end
    s = [num2str(k/length(flist)*100,'%.0f'),'% surface altitude data loaded, thanks for your patience'];
    waitbar(k/length(flist),w_l,s);
end
close(w_l);clear w_l;
end