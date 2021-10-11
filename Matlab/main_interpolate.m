% main script to interpolate data.

addpath('../Data/raw_data/');
load('fr_800_with_coordinate_interpolated.mat');

size(interpolate_brain_image)
interpolate_brain_image2 = interpolateDeadElectrodes(interpolate_brain_image);
save('../Data/raw_data/fr_800_with_coordinate_interpolated2.mat', 'interpolate_brain_image2');
