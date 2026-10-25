function [smb_out, x_out, y_out] = loadSMBGreenland(smb_root, model_request)
%LOADSMB This function loads various SMB products from publicly available
%data set, and output to standardized multi-decadal mean format
%
%   Input:
%       smb_root [string]: 
%           root directory to the SMB reanalysis models
%       model_request [string]: requested regional climate model
%           currently supporting:
%               "RACMO"
%               "HIRHAM"
%               "MAR"
%   Output:
%       smb_out  [double]: MxN SMB in, meter water equiv. per year
%       x_out[double]: MxN grid of x coordinates in 3031 projection
%       y_out[double]: MxN grid of y coordinates in 3031 projection
%
    root_dir = smb_root;
    rho_w = 1e3; % kg m^-3
    switch model_request
        case "RACMO"
            smb_path = 'RACMO/GRN11_RACMO24p1_data/smbgl_monthlyS_GRN11_RACMO24p1_200601_201512.nc';
            smb_path = [root_dir, smb_path];
            smb  = double(ncread(smb_path,'smbgl')); % unit: kg m-2
            day = double(ncread(smb_path,'time')); % days since 1950-01-01 00:00:00.0
            lat  = double(ncread(smb_path,'lat'));
            lon  = double(ncread(smb_path,'lon'));

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
            smb_filtered  = smb(:,:,1,keep_idx);

            smb_year_sum = cell(1, numel(complete_years));
            % get the sum for each year
            for ii = 1:numel(complete_years)
                year_idx = year(time_filtered)==complete_years(ii);
                smb_year_sum{ii} = sum(smb_filtered(:,:,1,year_idx),4);
            end
            smb_annual = cat(3, smb_year_sum{:});
            smb_annual = smb_annual./rho_w; % convert to m.eq/a
            
            % flatten to output
            smb_out   = smb_annual; 
            x_out = x;
            y_out = y;


        case "HIRHAM"
            smb_path  = 'HIRHAM/GREENLAND/ERAI/ANNUAL/DMI-HIRHAM5_GL2_ERAI_1980_2014_SMB_YM.nc';
            smb_path = [root_dir, smb_path];

            smb = ncread(smb_path, 'smb');
            smb = squeeze(smb);

            % load the coordinate
            lat = double(ncread(smb_path,'lat'));
            lon = double(ncread(smb_path,'lon'));
            original_shape = size(lat);
            [x_flat,y_flat] = ll2psn(lat(:), lon(:));
            x = reshape(x_flat, original_shape);
            y = reshape(y_flat, original_shape);

            % output
            % convert from mm water equivalent to meters
            smb_out = smb/1e3;
            x_out = x;
            y_out = y;


        case "MAR"
            % this is a daily dataset
            smb_path = 'MAR/data/MARv3_11_20km_daily_ERA_10km_1979_2019.nc';
            smb_path = [root_dir, smb_path];

            % load
            day = double(ncread(smb_path,'TIME'));

            % using lat, lon instead of the X10_153 and Y21_288 variables
            % so no need to transpose
            %x = double(ncread(smb_path,'X10_153'));
            %y = double(ncread(smb_path,'Y21_288'));
            lat = double(ncread(smb_path,'LAT'));
            lon = double(ncread(smb_path,'LON'));
            original_shape = size(lat);
            [x_flat,y_flat] = ll2psn(lat(:), lon(:));
            x = reshape(x_flat, original_shape);
            y = reshape(y_flat, original_shape);

            % notice that SMB unit in the NetCDF is incorrect!
            % it is mm. w.e./day
            smb  = double(ncread(smb_path,'SMB')); % mm. w.e./day

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
            smb_filtered  = smb(:,:,1,keep_idx);

            smb_year_sum = cell(1, numel(complete_years));
            % get the sum for each year
            for ii = 1:numel(complete_years)
                % find all the monthly data of a given year
                year_idx = year(time_filtered)==complete_years(ii);
                smb_year_sum{ii} = sum(smb_filtered(:,:,1,year_idx), 4, 'omitnan');
            end
            smb_annual_mean = cat(3, smb_year_sum{:});

            % convert to m.w.e/yr
            smb_annual_mean = smb_annual_mean/1e3;

            %[x_out, y_out] = meshgrid(x,y);

            % output
            smb_out   = smb_annual_mean;
            x_out = x;
            y_out = y;


        otherwise
            error("Unknown model specified!")
    end

end

