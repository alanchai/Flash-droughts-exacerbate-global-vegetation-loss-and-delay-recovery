clc;
clear;

% Load data files
% load('global_NDVI_mean.mat');
% load('global_NDVI_std.mat');
% load('global_NDVI_1982_2022.mat');
% load('Drought_onset_start.mat');
% load('Drought_recovery_end.mat');


[lat_size, lon_size, num_time_points] = size(Drought_recovery_end);
lag_pre = zeros(lat_size, lon_size, num_time_points);
lag_post = zeros(lat_size, lon_size, num_time_points);
Drought_GPP_with_lag = zeros(lat_size, lon_size, num_time_points); 
Drought_GPP_percentage_with_lag = zeros(lat_size, lon_size, num_time_points);
GPP_recovery_time = zeros(lat_size, lon_size, num_time_points);
Drought_GPP_no_lag = zeros(lat_size, lon_size, num_time_points); 
Drought_GPP_percentage_no_lag = zeros(lat_size, lon_size, num_time_points); 
GPP_before_drought_with_lag = zeros(lat_size, lon_size, num_time_points);
GPP_after_drought_with_lag = zeros(lat_size, lon_size, num_time_points);
GPP_before_drought_no_lag = zeros(lat_size, lon_size, num_time_points);
GPP_after_drought_no_lag = zeros(lat_size, lon_size, num_time_points);

% Load land mask data
load('land.mat');
land_data = land;

[valid_lat_indices, valid_lon_indices] = find(land_data ==1);
num_valid_points = length(valid_lat_indices);
for point_idx = 1:num_valid_points % Process all valid points
    lon_idx = valid_lon_indices(point_idx);
    lat_idx = valid_lat_indices(point_idx);
    onset_dates = Drought_onset_start(lon_idx, lat_idx, :);
    num_onset_dates = size(onset_dates, 3);
    onset_dates_vector = squeeze(onset_dates(:));
    [valid_onset_dates, valid_onset_indices] = find(onset_dates_vector > 0);
    num_valid_onset_dates = length(valid_onset_dates);
    for valid_idx = 1:num_valid_onset_dates
        onset_start = valid_onset_dates(valid_idx);
        onset_end = Drought_recovery_end(lon_idx, lat_idx, onset_start);
        NDVI_series = lastt(lon_idx, lat_idx, onset_start:onset_end);
        % Spline fitting for NDVI data
        start_cycle = onset_start / 73;
        start_cycle_idx = fix(start_cycle);
        start_day_offset = onset_start - start_cycle_idx * 73;
        end_cycle = onset_end / 73;
        end_cycle_idx = fix(end_cycle);
        end_day_offset = onset_end - end_cycle_idx * 73;
        cycle_start_adjusted = start_cycle_idx + 1;
        num_days = end_day_offset - start_day_offset + 1;
        NDVI_data_matrix = zeros(num_days, 41);
        for year_idx = 1:41
            yearly_start = (year_idx - 1) * 73 + 1;
            yearly_end = yearly_start + 72;
            yearly_data = lastt(lon_idx, lat_idx, yearly_start:yearly_end);
            NDVI_data_matrix(:, year_idx) = squeeze(yearly_data(:));
        end
        NDVI_adjusted = zeros(num_days, 1);
        for day = start_day_offset:end_day_offset
            y_data = NDVI_data_matrix(day, :)';
            x_data = (1:41)';
            valid_data_indices = ~isnan(y_data);
            if any(valid_data_indices)
                y_data = y_data(valid_data_indices);
                x_data = x_data(valid_data_indices);
                spline_fit = fit(x_data, y_data, 'smoothingspline', 'SmoothingParam', 0.1);
                NDVI_adjusted(day - start_day_offset + 1) = NDVI_series(day) - spline_fit(cycle_start_adjusted);
            else
                NDVI_adjusted(day - start_day_offset + 1) = NaN;
            end
        end
        count = sum(NDVI_adjusted <= 0);
        lag_pre(lon_idx, lat_idx, onset_start) = onset_end - count + 1;
        lag_end_adjusted = onset_end + 1;
        NDVI_next_series = lastt(lon_idx, lat_idx, lag_end_adjusted:2993);
        NDVI_previous_series = lastt(lon_idx, lat_idx, onset_end:2992);
        NDVI_change = NDVI_next_series - NDVI_previous_series;
        count = 0;
        for i = 1:length(NDVI_change)
            if NDVI_change(i) <= 0
                count = count + 1;
            else
                break;
            end
        end
        lag_post(lon_idx, lat_idx, onset_start) = onset_end + count;
        NDVI_data_yearly = zeros(num_days, 41);
        for year_idx = 1:41
            yearly_start = (year_idx - 1) * 73 + 1;
            yearly_end = yearly_start + 72;
            yearly_data = lastt(lon_idx, lat_idx, yearly_start:yearly_end);
            NDVI_data_yearly(:, year_idx) = squeeze(yearly_data(:));
        end
        NDVI_fit_adjusted = zeros(73, 1);
        for day_idx = 1:73
            y_data = NDVI_data_yearly(day_idx, :)';
            x_data = (1:41)';
            valid_data_indices = ~isnan(y_data);
            if any(valid_data_indices)
                y_data = y_data(valid_data_indices);
                x_data = x_data(valid_data_indices);
                spline_fit = fit(x_data, y_data, 'smoothingspline', 'SmoothingParam', 0.1);
                NDVI_fit_adjusted(day_idx) = spline_fit(cycle_start_adjusted);
            else
                NDVI_fit_adjusted(day_idx) = NaN;
            end
        end
        recovery_start = lag_post(lon_idx, lat_idx, onset_start) + 1;
        if recovery_start >= 2993%total number of pentads
            GPP_recovery_time(lon_idx, lat_idx, onset_start) = NaN;
        else
            NDVI_difference_series = zeros(2993 - recovery_start + 1, 1);
            for itr = recovery_start:2993
                day_offset = itr - recovery_start + 1;
                NDVI_difference_series(day_offset) = lastt(lon_idx, lat_idx, itr) - NDVI_fit_adjusted(mod(itr, 73) + 1);
            end
            index = find(all(NDVI_difference_series(1:end-2) >= 0));
            if isempty(index)
                GPP_recovery_time(lon_idx, lat_idx, onset_start) = 2993;
            else
                GPP_recovery_time(lon_idx, lat_idx, onset_start) = lag_post(lon_idx, lat_idx, onset_start) + index(1);
            end
        end
        historical_GPP = nanmean(NDVI_series);
        GPP_after_drought_no_lag(lon_idx, lat_idx, onset_start) = historical_GPP;
        GPP_after_drought_with_lag(lon_idx, lat_idx, onset_start) = GPP_after_FDD; 
        Drought_GPP_with_lag(lon_idx, lat_idx, onset_start) = Drought_effect;
        Drought_GPP_percentage_with_lag(lon_idx, lat_idx, onset_start) = Drought_percentage;
    end
end

save('Cai_Drought_on_GPP_with_lag.mat', 'lag_pre', 'lag_post', 'Drought_GPP_with_lag', 'Drought_GPP_percentage_with_lag', 'GPP_recovery_time', 'Drought_GPP_no_lag', 'Drought_GPP_percentage_no_lag', 'GPP_before_drought_with_lag', 'GPP_after_drought_with_lag');