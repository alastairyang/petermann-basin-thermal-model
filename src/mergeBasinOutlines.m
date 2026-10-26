function [xi, yi, basins] = mergeBasinOutlines(filename, targetIDs)
%MERGEBASINOUTLINES Extract and merge basin outlines without downsampling.
%
% Usage:
%   [xi, yi, basins] = mergeBasinOutlines( ...
%       'grndrainagesystems_ekholm.txt', ["1.1", "1.2"]);
%
% Inputs:
%   filename  - Input text file containing an END OF HEADER marker.
%   targetIDs - Numeric IDs, string IDs, or a cell array of character IDs.
%
% Outputs:
%   xi, yi    - Closed, ordered outline of the merged basins, with all
%               shared input vertices removed.
%   basins    - Structure array containing the original projected outlines:
%                 basins(k).id
%                 basins(k).xi
%                 basins(k).yi
%
% Assumptions:
%   - Each basin's records describe its boundary in contour order.
%   - Shared vertices have exactly matching projected coordinates.
%   - The merged basins form one connected region without holes.
%
% Dependencies:
%   - ll2psn, using the same projection settings as your original script.
%   - MATLAB polyshape.
%
% No downsampling is performed.

    % Validate and normalize IDs

    if ~isnumeric(targetIDs)
        targetIDs = str2double(string(targetIDs));
    end

    validateattributes(targetIDs, {'numeric'}, ...
        {'vector', 'nonempty', 'real', 'finite'}, ...
        mfilename, 'targetIDs');

    targetIDs = unique(double(targetIDs(:)), 'stable');
    nBasins = numel(targetIDs);

    % Open file and skip header

    fid = fopen(filename, 'r');

    if fid == -1
        error('mergeBasinOutlines:OpenFailed', ...
            'Could not open file: %s', filename);
    end

    cleanup = onCleanup(@() fclose(fid));

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
        error('mergeBasinOutlines:MissingHeader', ...
            'The file does not contain an END OF HEADER marker.');
    end

    % Read subsystem ID, latitude, and longitude

    data = textscan(fid, '%f %f %f', ...
        'MultipleDelimsAsOne', true);

    ids = data{1};
    lat = data{2};
    lon = data{3};

    if numel(ids) ~= numel(lat) || numel(ids) ~= numel(lon)
        error('mergeBasinOutlines:InvalidData', ...
            'The input columns have inconsistent lengths.');
    end

    % Extract individual projected outlines

    basins = repmat(struct('id', [], 'xi', [], 'yi', []), ...
        nBasins, 1);

    points = cell(nBasins, 1);
    shared = cell(nBasins, 1);

    for k = 1:nBasins
        mask = abs(ids - targetIDs(k)) < 1e-8;

        if ~any(mask)
            error('mergeBasinOutlines:MissingBasin', ...
                'No records found for basin ID %.15g.', targetIDs(k));
        end

        [xk, yk] = ll2psn(lat(mask), lon(mask));

        xk = xk(:);
        yk = yk(:);

        % Preserve the original projected coordinates.
        basins(k).id = targetIDs(k);
        basins(k).xi = xk;
        basins(k).yi = yk;

        pk = [xk, yk];

        if any(~isfinite(pk(:)))
            error('mergeBasinOutlines:InvalidCoordinates', ...
                'Basin %.15g contains nonfinite coordinates.', ...
                targetIDs(k));
        end

        if size(unique(pk, 'rows'), 1) < 3
            error('mergeBasinOutlines:InvalidOutline', ...
                'Basin %.15g has fewer than three distinct vertices.', ...
                targetIDs(k));
        end

        points{k} = pk;
        shared{k} = false(size(pk, 1), 1);
    end

    % Identify vertices shared by different basins
    % Repeated vertices within a single basin, including its closing
    % vertex, do not count as shared.

    for k = 1:nBasins-1
        for j = k+1:nBasins
            shared{k} = shared{k} | ...
                ismember(points{k}, points{j}, 'rows');

            shared{j} = shared{j} | ...
                ismember(points{j}, points{k}, 'rows');
        end
    end

    sharedPoints = zeros(0, 2);

    for k = 1:nBasins
        sharedPoints = [sharedPoints; ...
            points{k}(shared{k}, :)]; 
    end

    sharedPoints = unique(sharedPoints, 'rows');

    % Merge polygon interiors to obtain the ordered outer boundary
    % Merge BEFORE deleting shared vertices: deleting them first would
    % destroy the topology needed to connect the remaining boundary arcs.

    merged = polyshape(points{1}(:, 1), points{1}(:, 2), ...
        'KeepCollinearPoints', true);

    for k = 2:nBasins
        current = polyshape(points{k}(:, 1), points{k}(:, 2), ...
            'KeepCollinearPoints', true);

        merged = union(merged, current, ...
            'KeepCollinearPoints', true);
    end

parts = regions(merged);

if isempty(parts)
    error('mergeBasinOutlines:EmptyUnion', ...
        'The polygon union is empty.');
end

partAreas = arrayfun(@(p) area(p), parts);
[largestArea, mainIndex] = max(partAreas);

otherAreas = partAreas;
otherAreas(mainIndex) = [];

% Only permit negligible disconnected fragments.
% This is an AREA tolerance, not a coordinate-matching tolerance.
relativeAreaTolerance = 1e-10;

if any(otherAreas > relativeAreaTolerance * largestArea)
    error('mergeBasinOutlines:DisconnectedBasins', ...
        ['The union contains %d disconnected regions. ', ...
         'Largest area: %.6g; total area of other regions: %.6g. ', ...
         'Check shared-boundary coordinates before merging.'], ...
        numel(parts), largestArea, sum(otherAreas));
end

if numel(parts) > 1
    warning('mergeBasinOutlines:SmallFragments', ...
        'Discarding %d negligible disconnected fragment(s).', ...
        numel(parts) - 1);
end

outerPolygon = parts(mainIndex);

% Extract only the exterior boundary
% Inner loops are not part of the requested exterior outline.

if outerPolygon.NumHoles > 0
    warning('mergeBasinOutlines:InteriorHoles', ...
        ['The union contains %d hole(s). Returning only its ', ...
         'exterior outline; inspect the holes for possible ', ...
         'shared-boundary mismatches.'], ...
        outerPolygon.NumHoles);
end

outerPolygon = rmholes(outerPolygon);

[xBoundary, yBoundary] = boundary(outerPolygon);
outline = [xBoundary(:), yBoundary(:)];

% A single exterior loop should contain no NaN separators.
if any(~isfinite(outline(:)))
    error('mergeBasinOutlines:UnexpectedBoundary', ...
        'The exterior boundary still contains multiple loops.');
end

% Remove the closing duplicate before filtering.
if size(outline, 1) > 1 && ...
        isequal(outline(1, :), outline(end, :))
    outline(end, :) = [];
end

% Remove shared input vertices
% This follows your original requirement to remove every shared point.

if ~isempty(sharedPoints)
    outline(ismember(outline, sharedPoints, 'rows'), :) = [];
end

% Remove consecutive duplicate vertices, if present.
if ~isempty(outline)
    keep = [true; any(diff(outline, 1, 1) ~= 0, 2)];
    outline = outline(keep, :);
end

if size(unique(outline, 'rows'), 1) < 3
    error('mergeBasinOutlines:TooFewVertices', ...
        'Fewer than three distinct exterior vertices remain.');
end

% Close and return the contour

outline(end+1, :) = outline(1, :);

xi = outline(:, 1);
yi = outline(:, 2);
end
