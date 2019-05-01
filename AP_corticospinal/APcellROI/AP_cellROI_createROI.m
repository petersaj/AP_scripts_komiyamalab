% draw ROI polygon
[cellMask,polygon.ROI{n_polygon}(:,1),polygon.ROI{n_polygon}(:,2)] = roipoly;

% make the polygon closed if it's not already
if polygon.ROI{n_polygon}(1,1) ~= polygon.ROI{n_polygon}(end,1)
    polygon.ROI{n_polygon} = [polygon.ROI{n_polygon} polygon.ROI{n_polygon}(1,:)];
end

% display just an outline when not selected
polygon.handles{n_polygon} = line(polygon.ROI{n_polygon}(:,1), polygon.ROI{n_polygon}(:,2));
set(polygon.handles{n_polygon}, 'ButtonDownFcn', {@selectROI, handles});

