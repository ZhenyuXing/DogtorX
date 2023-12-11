function [air] = load_air(fname, minlat, maxlat, minlon, maxlon)
% LOAD AIR: This function is to derive column of dry air data from TROPOMI file which has been labeled with suspect plumes.
%   Input argugments should include the .nc file name and [minlat maxlat][minlon maxlon].
%   Output argument is the column of dry air in the region of interest on the corresponding day.
[grid_lon, grid_lat] = meshgrid(minlon:0.05:maxlon,minlat:0.05:maxlat);
[a,b] = size(grid_lat);
air = zeros(a,b);
xair = ncread(fname,'/PRODUCT/SUPPORT_DATA/INPUT_DATA/dry_air_subcolumns');
xair = mean(xair,1,'omitnan');
xair = squeeze(xair);
latitude = ncread(fname,'/PRODUCT/latitude');
longitude = ncread(fname,'/PRODUCT/longitude');
w_l = waitbar(0,'Extracting air data, please wait...');
for j = 1:b
    for i = 1:a
        poly_lat = [minlat+0.05*(i-1) minlat+0.05*(i-1) minlat+0.05*i minlat+0.05*i minlat+0.05*(i-1)];
        poly_lon = [minlon+0.05*(j-1) minlon+0.05*j minlon+0.05*j minlon+0.05*(j-1) minlon+0.05*(j-1)];
        in = inpolygon(longitude, latitude, poly_lon, poly_lat);
        air_raw = xair(in);
        air(i,j) = mean(air_raw(:),'omitnan');
    end
s = [num2str(j/b*100,'%.0f'),'% air data loaded, thanks for your patience'];
waitbar(j/b,w_l,s);    
end
close(w_l);clear w_l;
end