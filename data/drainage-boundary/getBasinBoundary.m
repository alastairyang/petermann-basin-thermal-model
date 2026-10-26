
% Merge 1.1 and 1.2, then downsample
[xi, yi, basins] = mergeBasinOutlines( ...
    'grndrainagesystems_ekholm.txt', ["1.1", "1.2"]);

% Subsample both groups
skip_coast    = 2;
skip_interior = 150;

% Subsample interior points and retain all coastal points
coast = expread('1.1-1.2-coast.exp');

% Column vectors for consistent indexing
xi = xi(:);
yi = yi(:);

in_coast = inpolygon(xi, yi, coast.x, coast.y);
interior = ~in_coast;

% Original indices of each group
idx_coast = find(in_coast);
idx_interior = find(interior);

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
coordToExp([xi_merged, yi_merged], '1.1-1.2-downsample-boundary.exp')