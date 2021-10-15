 %clear
%% HCP rfMRI
% addpath('D:\spatiotemporal patterns\HCP');
% load('aal2_2mm_mask.mat');
% mask = aal2_2mm_mask;
% 
% addpath('D:\spatiotemporal patterns\HCP\rfMRI\101006\MNINonLinear\Results\rfMRI_REST1_LR');
% fMRI = load_nii('rfMRI_REST1_LR.nii.gz');
% BrainImg = fMRI.img;
% 
% addpath('D:\spatiotemporal patterns\HCP\rfMRI\101006');
% load('VelField_3D_10_100_1_constrained.mat');

%% HCP tfMRI
% addpath('D:\spatiotemporal patterns\HCP');
% load('aal2_2mm_mask.mat');
% mask = aal2_2mm_mask;
% addpath('D:\spatiotemporal patterns\HCP\tfMRI_EMOTION\100206\MNINonLinear\Results\tfMRI_EMOTION_LR');
% fMRI = load_nii('tfMRI_EMOTION_LR.nii.gz');
% BrainImg = fMRI.img;
% %BrainImg = (BrainImg - mean(BrainImg,4))./std(BrainImg,1,4);
% addpath('D:\spatiotemporal patterns\HCP\tfMRI_EMOTION\100206');
% %load('VelField_3D_176_HS.mat');
% load('VelField_3D_176_10_10_constrained.mat');

%% UK biobank rfMRI
% load('D:\spatiotemporal patterns\UK biobank\aal2_mask.mat');
% mask = aal2_mask;
% 
% dir ='D:\spatiotemporal patterns\UK biobank\UKB_1\2333381\';
% fMRI = load_nii([dir,'sFunImg_3mmStdSpace.nii.gz']);
% BrainImg = fMRI.img;
% %BrainImg = BrainImg*500 +10000;
% addpath('D:\spatiotemporal patterns\UK biobank\UKB_1\2333381');
% %load('VelField_3D_340_1000_HS.mat');%VelField = VelField_3d;
% load([dir,'VelField_3D_490_100_1_constrained2.mat']);

[M,N,S,T] = size(Ux);
% %Make a mask with size of M*N*S
% C_mask = zeros([M N S]);
% for i = 1:M
%     for j = 1:N
%         for s = 1:S
%             if mask(i,j,s)==1 && mask(i,j+1,s)==1 && mask(i+1,j,s)==1 && mask(i+1,j+1,s)==1....
%                     && mask(i,j,s+1)==1 && mask(i,j+1,s+1)==1 && mask(i+1,j,s+1)==1 && mask(i+1,j+1,s+1)==1
%                 C_mask(i,j,s) = 1;
%             end
%         end
%     end
% end
BrainImg = interpolate_brain_image;
BrainImg(BrainImg==0)=nan;
[X,Y] = meshgrid(1:N,1:M);
X=X+0.5; Y=Y+0.5;
%T = 9;
z = 30;
H(T) = struct('cdata',[],'colormap',[]);dd=1;
for t = 1:T
    imagesc(BrainImg(:,:,z,t));hold on % .*mask(:,:,z)
    quiver(X(1:dd:end,1:dd:end),Y(1:dd:end,1:dd:end),...
        Ux(1:dd:end,1:dd:end,z,t),...
        Uy(1:dd:end,1:dd:end,z,t),3); hold off
    axis([0 N 0 M]);
    pause(0.5)
    %quiver(X,Y,VelField(:,:,z,1,t),VelField(:,:,z,2,t),5); hold off
    drawnow
    H(t) = getframe;
end
%[M,N,O,~,T] = size(VelField);
% [X,Y,Z] = meshgrid(1:N,1:M,1:S);
% X=X+0.5; Y=Y+0.5;Z=Z+0.5;
% %T = 9;
% z = 30;
% H(T) = struct('cdata',[],'colormap',[]);dd=5;
% for t = 1:T
%     %imagesc(BrainImg(:,:,z,t).*mask(:,:,z));hold on % .*mask(:,:,z)
%     quiver3(X(1:dd:end,1:dd:end,1:dd:end),Y(1:dd:end,1:dd:end,1:dd:end),Z(1:dd:end,1:dd:end,1:dd:end),...
%         Ux(1:dd:end,1:dd:end,1:dd:end,t).*C_mask(1:dd:end,1:dd:end,1:dd:end),...
%         Uy(1:dd:end,1:dd:end,1:dd:end,t).*C_mask(1:dd:end,1:dd:end,1:dd:end),...
%         Uz(1:dd:end,1:dd:end,1:dd:end,t).*C_mask(1:dd:end,1:dd:end,1:dd:end),1); %hold off
%     axis([0 N 0 M 0 S]);
%     %quiver(X,Y,VelField(:,:,z,1,t),VelField(:,:,z,2,t),5); hold off
%     drawnow
%     H(t) = getframe;
% end
% %axis([0 N 0 M]);hold on
% % movie(H,1,0.5);