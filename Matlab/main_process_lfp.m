%% main script to process data and compute velocity.
% including : extract phase, rexamine phase, interpolate, compute velocity.
addpath('../Data/raw_data/');
load('lfp_200hz.mat');
disp('size:'); 
brain_image = permute(brain_image, [2, 1, 3, 4]); % (Y, X, Z, T)
size(brain_image)
%% reduce resolution
disp('reduce resolution...');
brain_image = reduce_resolution(brain_image);
disp('After reduce resolution, size')
size(brain_image)

% some channels have abnormal LFPs.
resample = true; 

% Find any channels with NaN or zero values
nanChans = any(isnan(brain_image(:,:,:, :)),4);
zeroChans = all(brain_image(:,:,:, :)==0, 4);
badChannels = find(nanChans | zeroChans);

if resample
    % Also exclude channel with abnormal LFPs
    abnormalChannels = find(any((brain_image<-65), 4));
    badChannels = union(abnormalChannels, badChannels);
end
%% Filter by band-pass filters, TODO, there are some bug.
% fLow = 8;
% fHigh = 12;
% Fs = 1000;
% tic
% brain_image = filterSignal(brain_image, fLow, fHigh, Fs);
% toc 
% ind = find(all(~isnan(brain_image), 4));
% ind = datasample(ind, 10);
% sz = size(brain_image, [1, 2, 3]);
% [x, y, z] = ind2sub(sz, ind);
% 
% figure;
% for i=1:10
%     subplot(10,1,i);
%     data = squeeze(brain_image(x(i), y(i), z(i), :));
%     plot(data);
% end
% saveas(gcf, 'figure/bandfilter_lfp.jpg');
% close;
%% find mask in HPC data;
% BrainImg = permute(BrainImg, [2, 1, 3, 4]);
% BrainImg(BrainImg==0)=nan;
% indB = find(all(~isnan(BrainImg), 4));
% szB = size(BrainImg, [1, 2, 3]);
% [yy, xx, zz] = ind2sub(szB, indB);
% 
% indb = find(all(~isnan(brain_image), 4));
% szb = size(brain_image, [1, 2, 3]);
% [y, x, z] = ind2sub(szb, indb);
% yratio = (max(yy)-min(yy)) / (max(y)-min(y));
% xratio = (max(xx)-min(xx)) / (max(x)-min(x));
% zratio = (max(zz)-min(zz)) / (max(z)-min(z));
% yy = round(yy / yratio);
% xx = round(xx / xratio);
% zz = round(zz / zratio);
% coords = [xx, yy, zz];
% coords = unique(coords, 'rows');
% %TODO: TO FINISH IT

%% Filter LFPs to a fixed frequedny， seems it is not a good choice.
% disp('Filtering waveforms...'); 
% tic
% timeDim = 4;
% Fs = 1000;
% % morleCfreqs = linspace(1, 20, 100);
% morleCfreqs = 10;
% morletParam = 7;
% wvcfs = squeeze(morletWaveletTransform(brain_image, Fs, morleCfreqs, morletParam, timeDim));
% toc
%% Visualize the raw lfps and phase of 10hz frequency
% ind = find(all(~isnan(brain_image), 4));
% ind = datasample(ind, 10);
% sz = size(brain_image, [1, 2, 3]);
% [x, y, z] = ind2sub(sz, ind);
% 
% figure;
% for i=1:10
%     subplot(10,1,i);
%     data = squeeze(brain_image(x(i), y(i), z(i), :));
%     plot(data);
% end
% saveas(gcf, 'figure/raw_lfp.jpg');
% close;
% figure;
% for i=1:10
%     subplot(10,1,i);
%     data = squeeze(angle(wvcfs(x(i), y(i), z(i), :)));
%     plot(data);
% end
% saveas(gcf, 'figure/phase_lfp.jpg');
% close;

%% Test the maximum difference between successive time points
% useAmplitude = false;
% downsampleScale = 1;
% if useAmplitude
%     maxRange = max(abs(wvcfs),[],timeDim) - min(abs(wvcfs),[],timeDim);
%     allDiff = abs(diff(abs(wvcfs),1,timeDim)) ./ ...
%         repmat(maxRange,1,1,size(wvcfs,3)-1,1);
% else
%     allDiff = anglesubtract(angle(wvcfs(:,:,2:end,:)), ...
%         angle(wvcfs(:,:,1:end-1,:)), 1);
%     allDiff = allDiff / pi;
% end
% maxDiff = prctile(abs(allDiff(:)), 99);
% fprintf('99th percentile of the fractional change between time steps is %0.2f.\n', maxDiff)
% if maxDiff > 0.1
%     disp('Change is >10%, results may be affected by low sampling frequency.')
%     if downsampleScale > 1
%         disp('Consider reducing downsampleScale in setParams.m to increase Fs.')
%     end
% else
%     disp('Change is <10%, sampling frequency is sufficient for data.')
%     disp('Consider increasing downsampleScale in setParams.m to improve speed.')
% end
% 
% clearvars allDiff

%% 
disp('interpolate');
tic
interpolate_brain_image = interpolate(brain_image, badChannels);
nanChans2 = any(isnan(interpolate_brain_image(:,:,:, :)),4);
zeroChans2 = all(interpolate_brain_image(:,:,:, :)==0, 4);
badChannels2 = find(nanChans2 | zeroChans2);
disp('nan channels size')
size(badChannels2)
save('../Data/raw_data/interpolate_lfp.mat', 'interpolate_brain_image');
toc
