 function [Ux,Uy,Uz,res] = velocity_field_constrained_v2(BrainImg,para_model,para_alg)
%VELOCITY_FIELD_CONSTRAINED extracts the velocity field from brain images by 
%a tensor model solved by GADIM/MSS
% Input:
%   BrainImg: an (M+1)*(N+1)*(S+1)*(T+1) tensor (BOLD-fMRI signal)
%             y-, x-, z-, t-direction
%   para_model: model parameters including ...
%       rho: the penalty parameter for the smooth regularization
%       tau: the penalty parameter for the time-continuity
%   para_alg: algorithm parameters including ...
%       alpha,beta: stepsizes
%       omega: relaxation
%       tol: tolerance for stopping rule
%       maxit: maximum number of iterations
% Output:
%   Ux,Uy,Uz: M*N*S*T tensors of the x-, y- and z-components of the velocity field
%   res: residuals
%   by Weiyang Ding @Fudan August 30, 2021
%   modified @Fudan October 13, 2021

% Estimate the partial derivatives
conv_ker = repmat([1,-1]./8,[2,1,2,2]);
Dx = convn(BrainImg,conv_ker,'valid');% d/dx
Dy = convn(BrainImg,permute(conv_ker,[2,1,3,4]),'valid');% d/dy
Dz = convn(BrainImg,permute(conv_ker,[1,3,2,4]),'valid');% d/dz
Dt = convn(BrainImg,permute(conv_ker,[1,3,4,2]),'valid');% d/dt
% Initialization
res_primal = zeros(15,1);
res_dual = zeros(12,1);
res = zeros(2,para_alg.maxit);
[M,N,S,T] = size(Dx);% M - z-axis, N - x-axis, S - y-axis, T - t-axis
Px = zeros([M,N,S,T]); Py = zeros([M,N,S,T]); Pz = zeros([M,N,S,T]);
Qxx = zeros([M,N,S,T]); Qxy = zeros([M,N,S,T]); Qxz = zeros([M,N,S,T]);
Qyx = zeros([M,N,S,T]); Qyy = zeros([M,N,S,T]); Qyz = zeros([M,N,S,T]);
Qzx = zeros([M,N,S,T]); Qzy = zeros([M,N,S,T]); Qzz = zeros([M,N,S,T]);
Rx = zeros([M,N,S,T]); Ry = zeros([M,N,S,T]); Rz = zeros([M,N,S,T]);
Lxx = zeros([M,N,S,T]); Lxy = zeros([M,N,S,T]); Lxz = zeros([M,N,S,T]); Lxt = zeros([M,N,S,T]);
Lyx = zeros([M,N,S,T]); Lyy = zeros([M,N,S,T]); Lyz = zeros([M,N,S,T]); Lyt = zeros([M,N,S,T]);
Lzx = zeros([M,N,S,T]); Lzy = zeros([M,N,S,T]); Lzz = zeros([M,N,S,T]); Lzt = zeros([M,N,S,T]);
% Compute the FFT coefficients
cx = conj(fft([1;-1;zeros(N-2,1)])).';
cy = conj(fft([1;-1;zeros(M-2,1)]));
cz = conj(fft(permute([1;-1;zeros(S-2,1)],[2,3,1])));
ct = conj(fft(permute([1;-1;zeros(T-2,1)],[2,3,4,1])));
Cxyzt = ones([M,N,S,T]).*(1+para_alg.beta*para_alg.beta);
Cxyzt = bsxfun(@plus,Cxyzt,abs(cx).^2.*(para_alg.beta*para_alg.beta));
Cxyzt = bsxfun(@plus,Cxyzt,abs(cy).^2.*(para_alg.beta*para_alg.beta));
Cxyzt = bsxfun(@plus,Cxyzt,abs(cz).^2.*(para_alg.beta*para_alg.beta));
Cxyzt = bsxfun(@plus,Cxyzt,abs(ct).^2.*(para_alg.beta*para_alg.beta));
% Prepare the mask
mask = repmat(BrainImg(1:M,1:N,1:S,T)==0,1,1,1,T);
for iter = 1:para_alg.maxit
    % U-update
    Temp = (Dx.*Px + Dy.*Py + Dz.*Pz + Dt)./(Dx.*Dx + Dy.*Dy + Dz.*Dz);
    Temp = sign(Temp).*min(abs(Temp),para_alg.alpha);
    Temp(isnan(Temp)) = 0;
    Ux = Px - Dx.*Temp;
    Uy = Py - Dy.*Temp;
    Uz = Pz - Dz.*Temp;
    Ux(mask) = 0;
    Uy(mask) = 0;
    Uz(mask) = 0;
    % V-update
    Temp = (para_model.rho*para_alg.alpha)./...
        sqrt(Qxx.*Qxx+Qxy.*Qxy+Qxz.*Qxz...
        +Qyx.*Qyx+Qyy.*Qyy+Qyz.*Qyz...
        +Qzx.*Qzx+Qzy.*Qzy+Qzz.*Qzz);
    Temp = 1- min(Temp,1);
    Vxx = Qxx.*Temp; Vxy = Qxy.*Temp; Vxz = Qxz.*Temp;
    Vyx = Qyx.*Temp; Vyy = Qyy.*Temp; Vyz = Qyz.*Temp;
    Vzx = Qzx.*Temp; Vzy = Qzy.*Temp; Vzz = Qzz.*Temp;
    Temp = Ux - Ux(:,[2:N,1],:,:) - Vxx; res_dual(1) = norm(Temp(:));
    Temp = Ux - Ux([2:M,1],:,:,:) - Vxy; res_dual(2) = norm(Temp(:));
    Temp = Ux - Ux(:,:,[2:S,1],:) - Vxz; res_dual(3) = norm(Temp(:));
    Temp = Uy - Uy(:,[2:N,1],:,:) - Vyx; res_dual(4) = norm(Temp(:));
    Temp = Uy - Uy([2:M,1],:,:,:) - Vyy; res_dual(5) = norm(Temp(:));
    Temp = Uy - Uy(:,:,[2:S,1],:) - Vyz; res_dual(6) = norm(Temp(:));
    Temp = Uz - Uz(:,[2:N,1],:,:) - Vzx; res_dual(7) = norm(Temp(:));
    Temp = Uz - Uz([2:M,1],:,:,:) - Vzy; res_dual(8) = norm(Temp(:));
    Temp = Uz - Uz(:,:,[2:S,1],:) - Vzz; res_dual(9) = norm(Temp(:));
    % W-update
    Wx = Rx./(1+para_model.tau*para_alg.alpha);
    Wy = Ry./(1+para_model.tau*para_alg.alpha);
    Wz = Rz./(1+para_model.tau*para_alg.alpha);
    Temp = Ux - Ux(:,:,:,[2:T,1]) - Wx; res_dual(10) = norm(Temp(:));
    Temp = Uy - Uy(:,:,:,[2:T,1]) - Wy; res_dual(11) = norm(Temp(:));
    Temp = Uz - Uz(:,:,:,[2:T,1]) - Wz; res_dual(12) = norm(Temp(:));
    % lambda-update
    Temp = Ux.*2-Px;
    Lxx = Lxx - (Temp-Temp(:,[2:N,1],:,:)-Vxx.*2+Qxx).*para_alg.beta;
    Lxy = Lxy - (Temp-Temp([2:M,1],:,:,:)-Vxy.*2+Qxy).*para_alg.beta;
    Lxz = Lxz - (Temp-Temp(:,:,[2:S,1],:)-Vxz.*2+Qxz).*para_alg.beta;
    Lxt = Lxt - (Temp-Temp(:,:,:,[2:T,1])-Wx.*2+Rx).*para_alg.beta;
    etax = real(fftn(ifftn(Lxx - Lxx(:,[N,1:N-1],:,:)...
        + Lxy - Lxy([M,1:M-1],:,:,:) ...
        + Lxz - Lxz(:,:,[S,1:S-1],:) ...
        + Lxt - Lxt(:,:,:,[T,1:T-1]))./Cxyzt));
    Lxx = (Lxx - (etax-etax(:,[2:N,1],:,:)).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    Lxy = (Lxy - (etax-etax([2:M,1],:,:,:)).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    Lxz = (Lxz - (etax-etax(:,:,[2:S,1],:)).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    Lxt = (Lxt - (etax-etax(:,:,:,[2:T,1])).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    Temp = Uy.*2-Py;
    Lyx = Lyx - (Temp-Temp(:,[2:N,1],:,:)-Vyx.*2+Qyx).*para_alg.beta;
    Lyy = Lyy - (Temp-Temp([2:M,1],:,:,:)-Vyy.*2+Qyy).*para_alg.beta;
    Lyz = Lyz - (Temp-Temp(:,:,[2:S,1],:)-Vyz.*2+Qyz).*para_alg.beta;
    Lyt = Lyt - (Temp-Temp(:,:,:,[2:T,1])-Wy.*2+Ry).*para_alg.beta;
    etay = real(fftn(ifftn(Lyx - Lyx(:,[N,1:N-1],:,:)...
        + Lyy - Lyy([M,1:M-1],:,:,:) ...
        + Lyz - Lyz(:,:,[S,1:S-1],:) ...
        + Lyt - Lyt(:,:,:,[T,1:T-1]))./Cxyzt));
    Lyx = (Lyx - (etay-etay(:,[2:N,1],:,:)).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    Lyy = (Lyy - (etay-etay([2:M,1],:,:,:)).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    Lyz = (Lyz - (etay-etay(:,:,[2:S,1],:)).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    Lyt = (Lyt - (etay-etay(:,:,:,[2:T,1])).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    Temp = Uz.*2-Pz;
    Lzx = Lzx - (Temp-Temp(:,[2:N,1],:,:)-Vzx.*2+Qzx).*para_alg.beta;
    Lzy = Lzy - (Temp-Temp([2:M,1],:,:,:)-Vzy.*2+Qzy).*para_alg.beta;
    Lzz = Lzz - (Temp-Temp(:,:,[2:S,1],:)-Vzz.*2+Qzz).*para_alg.beta;
    Lzt = Lzt - (Temp-Temp(:,:,:,[2:T,1])-Wz.*2+Rz).*para_alg.beta;
    etaz = real(fftn(ifftn(Lzx - Lzx(:,[N,1:N-1],:,:)...
        + Lzy - Lzy([M,1:M-1],:,:,:) ...
        + Lzz - Lzz(:,:,[S,1:S-1],:) ...
        + Lzt - Lzt(:,:,:,[T,1:T-1]))./Cxyzt));
    Lzx = (Lzx - (etaz-etaz(:,[2:N,1],:,:)).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    Lzy = (Lzy - (etaz-etaz([2:M,1],:,:,:)).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    Lzz = (Lzz - (etaz-etaz(:,:,[2:S,1],:)).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    Lzt = (Lzt - (etaz-etaz(:,:,:,[2:T,1])).*(para_alg.beta*para_alg.beta))./(1+para_alg.beta*para_alg.beta);
    % P-update
    Temp = Ux-Px+etax.*para_alg.beta; res_primal(1) = norm(Temp(:));
    Px = Px + Temp.*para_alg.omega;
    Temp = Uy-Py+etay.*para_alg.beta; res_primal(2) = norm(Temp(:));
    Py = Py + Temp.*para_alg.omega;
    Temp = Uz-Pz+etaz.*para_alg.beta; res_primal(3) = norm(Temp(:));
    Pz = Pz + Temp.*para_alg.omega;
    % Q-update
    Temp = Vxx-Qxx-Lxx.*para_alg.beta; res_primal(4) = norm(Temp(:));
    Qxx = Qxx + Temp.*para_alg.omega;
    Temp = Vxy-Qxy-Lxy.*para_alg.beta; res_primal(5) = norm(Temp(:));
    Qxy = Qxy + Temp.*para_alg.omega;
    Temp = Vxz-Qxz-Lxz.*para_alg.beta; res_primal(6) = norm(Temp(:));
    Qxz = Qxz + Temp.*para_alg.omega;
    Temp = Vyx-Qyx-Lyx.*para_alg.beta; res_primal(7) = norm(Temp(:));
    Qyx = Qyx + Temp.*para_alg.omega;
    Temp = Vyy-Qyy-Lyy.*para_alg.beta; res_primal(8) = norm(Temp(:));
    Qyy = Qyy + Temp.*para_alg.omega;
    Temp = Vyz-Qyz-Lyz.*para_alg.beta; res_primal(9) = norm(Temp(:));
    Qyz = Qyz + Temp.*para_alg.omega;
    Temp = Vzx-Qzx-Lzx.*para_alg.beta; res_primal(10) = norm(Temp(:));
    Qzx = Qzx + Temp.*para_alg.omega;
    Temp = Vzy-Qzy-Lzy.*para_alg.beta; res_primal(11) = norm(Temp(:));
    Qzy = Qzy + Temp.*para_alg.omega;
    Temp = Vzz-Qzz-Lzz.*para_alg.beta; res_primal(12) = norm(Temp(:));
    Qzz = Qzz + Temp.*para_alg.omega;
    % R-update
    Temp = Wx-Rx-Lxt.*para_alg.beta; res_primal(13) = norm(Temp(:));
    Rx = Rx + Temp.*para_alg.omega;
    Temp = Wy-Ry-Lyt.*para_alg.beta; res_primal(14) = norm(Temp(:));
    Ry = Ry + Temp.*para_alg.omega;
    Temp = Wz-Rz-Lzt.*para_alg.beta; res_primal(15) = norm(Temp(:));
    Rz = Rz + Temp.*para_alg.omega;
    % Check convergence
    res(1,iter) = norm(res_primal,'inf');
    res(2,iter) = norm(res_dual,'inf');
    disp(['Step ',num2str(iter),': ',num2str(norm(res(:,iter),'inf'))]);
    if norm(res(:,iter),'inf') < para_alg.tol
        break
    end
end

end

