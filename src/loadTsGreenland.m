function [Ts_out, x_out, y_out] = loadTsGreenland(ts_root,model_request)
%LOADSMB This function loads various Ts products from publicly available
%data set, and output to standardized multi-decadal mean format
%
%   Input:
%       ts_root [string]: 
%           root directory to the surface temperature reanalysis models
%       model_request [string]: requested regional climate model
%           currently supporting:
%               "RACMO"
%               "HIRHAM"
%               "MAR"
%   Output:
%       Ts_out  [double]: MxN surface or near-surface temperature in Kelvin
%       x_out[double]: MxN grid of x coordinates in 3031 projection
%       y_out[double]: MxN grid of y coordinates in 3031 projection
%
    root_dir = ts_root;

    switch model_request
        case "RACMO"
            % RACMO only provides near surface temperature (tas)
            Ts_path = 'RACMO/GrIS_RACMO24p1_data/tas_monthlyA_GRN11_RACMO24p1_200601_201512.nc';
            Ts_path = [root_dir, Ts_path];
            Ts  = double(ncread(Ts_path,'tas')); % unit: kg m-2
            day = double(ncread(Ts_path,'time')); % days since 1950-01-01 00:00:00.0
            lat  = double(ncread(Ts_path,'lat'));
            lon  = double(ncread(Ts_path,'lon'));

            % note that the correct way for coord transformation is 
            % to flatten. Do not use the original matrix directly
            % even tho ll2psn can take in matrix
            original_shape = size(lat);
            [x_flat,y_flat] = ll2psn(lat(:), lon(:));
            x = reshape(x_flat, original_shape);
            y = reshape(y_flat, original_shape);

            % get multi-decadal mean 
            % by first breaking into years
            reference_day = datetime(1950,1,1);
            time = reference_day + days(day);

            % Extract years from the datetime array
            years = year(time);
            
            % Count instances for each year
            [unique_years, ~, idx] = unique(years);
            counts = accumarray(idx, 1);
            
            % since this is monthly dataset, we anticipate that every year
            % to have 12 instances. Filter out the once have <12 month
            complete_years = unique_years(counts == 12);

            % Filter the time array to keep only complete years
            keep_idx = ismember(year(time), complete_years);
            time_filtered = time(keep_idx);
            Ts_filtered  = Ts(:,:,1,keep_idx);

            Ts_year_sum = cell(1, numel(complete_years));
            % get the sum for each year
            for ii = 1:numel(complete_years)
                year_idx = year(time_filtered)==complete_years(ii);
                Ts_year_sum{ii} = sum(Ts_filtered(:,:,1,year_idx),4);
            end
            % this is a monthly dataset, divide 12 months
            % to get average temperature
            Ts_annual_mean = cat(3, Ts_year_sum{:}) / 12;
            
            % flatten to output
            Ts_out   = Ts_annual_mean; 
            x_out = x;
            y_out = y;


        case "HIRHAM"
            Ts_path  = 'HIRHAM/GREENLAND/ERAI/DAILY/DMI-HIRHAM5_GL2_ERAI_1980_2014_TSURF_DM.nc.gz';
            Ts_path = [root_dir, Ts_path];

            tsurf = ncread(Ts_path, 'tsurf');

            day = double(ncread(Ts_path,'time'));

            % load the coordinate
            lat = double(ncread(Ts_path,'lat'));
            lon = double(ncread(Ts_path,'lon'));
            original_shape = size(lat);
            [x_flat,y_flat] = ll2psn(lat(:), lon(:));
            x = reshape(x_flat, original_shape);
            y = reshape(y_flat, original_shape);

            % reference day
            reference_day = datetime(1980, 1, 1);
            time = reference_day + days(day);

            % Extract years from the datetime array
            years = year(time);
            
            % Count instances for each year
            [unique_years, ~, idx] = unique(years);
            counts = accumarray(idx, 1);
            
            % since this is daily dataset, we anticipate that every year
            % to have 365 or 366 instances. Filter out the once have >=365
            % days
            complete_years = unique_years(counts >= 365);

            % Filter the time array to keep only complete years
            keep_idx = ismember(year(time), complete_years);
            time_filtered = time(keep_idx);
            tsurf_filtered  = tsurf(:,:,1,keep_idx);

            tsurf_year_sum = cell(1, numel(complete_years));
            % get the sum for each year
            for ii = 1:numel(complete_years)
                % find all the monthly data of a given year
                year_idx = year(time_filtered)==complete_years(ii);
                tsurf_year_sum{ii} = sum(tsurf_filtered(:,:,1,year_idx), 4, 'omitnan');
            end
            tsurf_annual_mean = cat(3, tsurf_year_sum{:});
            % assume one year has 365 days for averaging
            tsurf_annual_mean = tsurf_annual_mean / 365;

            % output
            Ts_out = tsurf_annual_mean;
            x_out = x;
            y_out = y;


        case "MAR"
            % this is a daily dataset
            Ts_path = 'MAR/data/MARv3_11_20km_daily_ERA_10km_1979_2019.nc';
            Ts_path = [root_dir, Ts_path];

            % load
            day = double(ncread(Ts_path,'TIME'));

            % using lat, lon instead of the X10_153 and Y21_288 variables
            % so no need to transpose
            %x = double(ncread(smb_path,'X10_153'));
            %y = double(ncread(smb_path,'Y21_288'));
            lat = double(ncread(Ts_path,'LAT'));
            lon = double(ncread(Ts_path,'LON'));
            original_shape = size(lat);
            [x_flat,y_flat] = ll2psn(lat(:), lon(:));
            x = reshape(x_flat, original_shape);
            y = reshape(y_flat, original_shape);

            ST2 = ncread(Ts_path, 'ST2');

            % reference day
            reference_day = datetime(1977, 9, 1);
            time = reference_day + days(day);

            % Extract years from the datetime array
            years = year(time);
            
            % Count instances for each year
            [unique_years, ~, idx] = unique(years);
            counts = accumarray(idx, 1);
            
            % since this is daily dataset, we anticipate that every year
            % to have 365 or 366 instances. Filter out the once have >=365
            % days
            complete_years = unique_years(counts >= 365);

            % Filter the time array to keep only complete years
            keep_idx = ismember(year(time), complete_years);
            time_filtered = time(keep_idx);
            ST2_filtered  = ST2(:,:,1,keep_idx);

            ST2_year_sum = cell(1, numel(complete_years));
            % get the sum for each year
            for ii = 1:numel(complete_years)
                % find all the monthly data of a given year
                year_idx = year(time_filtered)==complete_years(ii);
                ST2_year_sum{ii} = sum(ST2_filtered(:,:,1,year_idx), 4, 'omitnan');
            end
            ST2_annual_mean = cat(3, ST2_year_sum{:});
            % assume one year has 365 days for averaging
            ST2_annual_mean = ST2_annual_mean / 365;
            % convert to Kelvin
            ST2_annual_mean = ST2_annual_mean + 273.15;

            %[x_out, y_out] = meshgrid(x,y);

            % output
            Ts_out   = ST2_annual_mean;
            x_out = x;
            y_out = y;


        otherwise
            error("Unknown model specified!")
    end

end

