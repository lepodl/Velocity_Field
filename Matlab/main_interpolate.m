% main script to interpolate data.

addpath('../Data/raw_data/');
load('lfp.mat');
load('DTI_voxel_network.mat');

size(brain_image)
brain_image(brain_image==0)=nan;
disp('interpolate once\n');
interpolate_brain_image = interpolateDeadElectrodes(brain_image, interp_idnex);
save('../Data/raw_data/lfp_interpolated.mat', 'interpolate_brain_image');
