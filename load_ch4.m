function [grid_lon,grid_lat,ch4,mean_ch4,coverage_count,flist_ex] = load_ch4(flist, minlat, maxlat, minlon, maxlon, qa_pass)
% LOAD CH4: This function is to derive dry air mole mixing ratio of methane from TROPOMI files.
%   Input argugments should include the .nc file names, [minlat maxlat][minlon maxlon], qa_value threshold.
%   Output argument are the XCH4 observations and coverage count of valid observations based on the regridded map in the region of interest.
[grid_lon, grid_lat] = meshgrid(minlon:0.05:maxlon,minlat:0.05:maxlat);
[a,b] = size(grid_lat);
ch4 = zeros(a,b,length(flist));
flist_ex = flist;
w_l = waitbar(0,'Loading CH_4 data, please wait...');
for k = 1:length(flist)
    try
        xch4 = ncread(flist(k).name,'/PRODUCT/methane_mixing_ratio_bias_corrected');
        latitude = ncread(flist(k).name,'/PRODUCT/latitude');
        longitude = ncread(flist(k).name,'/PRODUCT/longitude');
        qa_value = ncread(flist(k).name,'/PRODUCT/qa_value');
        roi_lat = [minlat minlat maxlat maxlat minlat];
        roi_lon = [minlon maxlon maxlon minlon minlon];
        [in_roi,on_roi] = inpolygon(longitude,latitude,roi_lon,roi_lat);
        in_roi = in_roi+on_roi;
        if any(in_roi(:))
            for j = 1:b
                parfor i = 1:a
                    poly_lat = [minlat+0.05*(i-1) minlat+0.05*(i-1) minlat+0.05*i minlat+0.05*i minlat+0.05*(i-1)];
                    poly_lon = [minlon+0.05*(j-1) minlon+0.05*j minlon+0.05*j minlon+0.05*(j-1) minlon+0.05*(j-1)];
                    in = inpolygon(longitude, latitude, poly_lon, poly_lat);
                    valid = (qa_value >= qa_pass);
                    ch4_raw = xch4(in & valid);
                    ch4(i,j,k) = mean(ch4_raw(:),'omitnan');
                end
            end
        else
            flist_ex(k).name = 'ROI not covered'; % Mark the file not covering the roi.
        end
    catch
        flist_ex(k).name = 'Void'; % Remove the file if data file is not readable.
        continue
    end
    s = [num2str(k/length(flist)*100,'%.0f'),'% CH_4 data loaded, thanks for your patience'];
    waitbar(k/length(flist),w_l,s);
end
close(w_l);clear w_l;
ch4(ch4 == 0) = nan; % replace filled 0 values with nan.
mean_ch4 = mean(ch4,3,"omitnan");
coverage_count = sum(~isnan(ch4),3);
coverage_count(coverage_count == 0) = nan;
end