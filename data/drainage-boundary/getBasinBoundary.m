% Input file
filename = 'grndrainagesystems_ekholm.txt';  % Replace with your file path

% Open the file
fid = fopen(filename, 'r');
if fid == -1
    error('Could not open file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

% Skip the header
headerFound = false;
while true
    line = fgetl(fid);

    if ~ischar(line)
        break
    end

    if strcmp(strtrim(line), 'END OF HEADER')
        headerFound = true;
        break
    end
end

if ~headerFound
    error('The file does not contain an END OF HEADER marker.');
end

% Read columns: subsystem ID, latitude, longitude
data = textscan(fid, '%f %f %f', ...
    'MultipleDelimsAsOne', true);

ids = data{1};
lat = data{2};
lon = data{3};

% Select records for ID 6.2
targetID = 6.2;
mask = abs(ids - targetID) < 1e-8;

lat_62 = lat(mask);
lon_62 = lon(mask);

% Combined array: column 1 = latitude, column 2 = longitude
latlon_62 = [lat_62, lon_62];

fprintf('Found %d records for ID %.1f.\n', nnz(mask), targetID);

% Save the extracted coordinates
results = table(lat_62, lon_62, ...
    'VariableNames', {'Latitude', 'Longitude'});
% writetable(results, 'drainage_6_2.csv');

% convert lat lon to EPSG 3413
[xi, yi] = ll2psn(lat_62, lon_62);

% Subsample interior points and retain all coastal points
coast = expread('6.2-coast.exp');

% Column vectors for consistent indexing
xi = xi(:);
yi = yi(:);

in_coast = inpolygon(xi, yi, coast.x, coast.y);
interior = ~in_coast;

% Original indices of each group
idx_coast = find(in_coast);
idx_interior = find(interior);

% Subsample both groups
skip_coast    = 3;
skip_interior = 200;

idx_coast_ds = idx_coast(1:skip_coast:end);
idx_interior_ds = idx_interior(1:skip_interior:end);

% Merge indices in their original contour order
idx_keep = sort([idx_coast_ds; idx_interior_ds]);

xi_merged = xi(idx_keep);
yi_merged = yi(idx_keep);

% Close the contour for the EXP convention, unless already closed
if ~isempty(xi_merged) && ...
        (xi_merged(end) ~= xi_merged(1) || yi_merged(end) ~= yi_merged(1))
    xi_merged(end+1, 1) = xi_merged(1);
    yi_merged(end+1, 1) = yi_merged(1);
end

% Plot the closed, connected contour
figure('Position', [2000, 300, 600, 600]);

plot(xi_merged, yi_merged, 'k.-', 'MarkerSize', 8);
axis equal;
grid on;
xlabel('x');
ylabel('y');
title('Subsampled coastal + interior contour');

% save to .exp file
coordToExp([xi_merged, yi_merged], '6.2-downsample-boundary.exp')