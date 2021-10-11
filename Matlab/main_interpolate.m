% main script to interpolate data.

addpath('../Data/raw_data/');
load('vmean_with_coordinate.mat');

size(brain_image)
interpolate_brain_image = interpolateDeadElectrodes(brain_image);
disp('interpolate once\n');
interpolate_brain_image = interpolateDeadElectrodes(interpolate_brain_image);
disp('interpolate twice\n');
save('../Data/raw_data/vmean_with_coordinate_interpolated2.mat', 'interpolate_brain_image');
