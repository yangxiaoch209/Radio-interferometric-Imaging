clear;
close all;

load('data_test_main_all.mat');

addpath data/
addpath data/img
addpath data/uvw
addpath data/vis
addpath lib/
addpath alg/

fprintf('Generating new data ... \n\n');
%% image loading
% [im, N, Ny, Nx] = util_read_image(image_file_name);

importdata W28.fits ;
im=double(ans);
im=im/(max(max(im)));
[Nx,Ny]=size(im);
im=imresize(im,1024./Nx);
[Nx,Ny]=size(im);N=Nx*Ny;

% figure; imagesc(real(log10((im))));set(gca,'YDir','normal');axis off;
% colorbar, axis image; caxis([-3.5 0]);%title('Recovered image');
% colormap(cubehelix);
% set(gca,'position',[-0.05,0.02,0.98,0.934])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

%% 参数
input_snr = 35 ; % noise level (on the measurements)

f_log_min=-3.5; 
f_log_max=0;
f_res=0.0005;
f_err=0.005;


param_sampling.p=0.3;
%% generate the sampling pattern
param_sampling.N = N; % number of pixels in the image
param_sampling.Nox = ox*Nx; % number of pixels in the image
param_sampling.Noy = oy*Ny; % number of pixels in the image

% [uw, vw, ~] = util_gen_sampling_pattern(sampling_pattern, param_sampling);
% 
% 
% if use_symmetric_fourier_sampling
%     uw = [uw; -uw];
%     vw = [vw; -vw];
% end

%%  SKA/VLA
load ('test_vla_10s.uvw.mat');
uw=uvw(:,1);vw=uvw(:,2);
% % MAX=max(max(real(uvw)));
% % uw=uw./MAX*pi;
% % vw=vw./MAX*pi;
% figure;
% plot(uw,vw,'.k');hold on;
% xlabel('u(λ)','FontName','Times New Roman','FontSize',16);
% ylabel('v(λ)','FontName','Times New Roman','FontSize',16);

%% compute weights
param_precond.N = N; % number of pixels in the image
param_precond.Nox = ox*Nx; % number of pixels in the image
param_precond.Noy = oy*Ny; % number of pixels in the image
[aWw] = util_gen_preconditioning_matrix(uw, vw, param_precond);

%% set the blocks structure
nWw = ones(length(uw), 1);
[u, v, ~, uvidx, aW, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);


%% measurement operator initialization

fprintf('Initializing the NUFFT operator\n\n');
tstart = tic;
[A, At, G, W, Gw] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
tend = toc(tstart);

fprintf('Initialization runtime: %ds\n\n', ceil(tend));
R = length(v);

y0 = cell(num_tests, 1);
y = cell(num_tests, 1);
aY = cell(num_tests, 1);
epsilon = cell(num_tests, 1);
epsilons = cell(num_tests, 1);
epsilonT = cell(num_tests, 1);
epsilonTs = cell(num_tests, 1);
y0f = cell(num_tests, 1);
yf = cell(num_tests, 1);
input_snr_v = cell(num_tests, 1);

%% generate noisy input data
for k = 1:num_tests
    [y0{k}, y{k}, y0f{k}, yf{k}, aY{k}, input_snr_v{k}] = util_gen_input_data(im, G, W, A, input_snr, use_different_per_block_input_snr, per_block_input_snr_delta, uvidx);
    
    if use_symmetric_fourier_sampling
        y0f{k} = [y0f{k}(uvidx{k}(1:end/2)); conj(y0f{k}(uvidx{k}(1:end/2)))];
        yf{k} = [yf{k}(uvidx{k}(1:end/2)); conj(yf{k}(uvidx{k}(1:end/2)))];
        for j = 1:R
            y{k}{j} = [y{k}{j}(uvidx{k}(1:end/2)); conj(y{k}{j}(uvidx{k}(1:end/2)))];
            y0{k}{j} = [y0{k}{j}(uvidx{k}(1:end/2)); conj(y0{k}{j}(uvidx{k}(1:end/2)))];
            aY{k}{j} = [aY{k}{j}(uvidx{k}(1:end/2)); conj(aY{k}{j}(uvidx{k}(1:end/2)))];
        end
    end
    
    
    [epsilonT{k}, epsilonTs{k}, epsilon{k}, epsilons{k}] = util_gen_L2_bounds(y{k}, ...
        input_snr, [], l2_ball_definition, stopping_criterion, use_same_stop_criterion, param_l2_ball);
end

%% gridding
if use_gridded_data
    [yT, T, Tw] = util_compute_gridding_data(y, G, Gw);
else
    T = G;
    Tw = Gw;
    yT = y;
end

%% sparsity operator definition

[Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
[Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);


%% compute the operator norm
fprintf('Computing operator norms ...\n');

fprintf('Natural W ...\n');
evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);

fprintf('No natural W ...\n');
Gw_ = spdiags(1./cell2mat(nW), 0, length(nWw), length(nWw)) * Gw;
Gwt_ = Gw_';
evl_no_natw = op_norm(@(x) Gw_ * A(x), @(x) At(Gwt_ * x), [Ny, Nx], 1e-6, 200, verbosity);
clear Gw_ Gwt_;

fprintf('Preconditioning ...\n');
evl_precond = op_norm(@(x) sqrt(cell2mat(aW)) .* (Gw * A(x)), @(x) At(Gw' * (sqrt(cell2mat(aW)) .* x)), [Ny, Nx], 1e-6, 200, verbosity);

evl_blocks = zeros(R, 1);

if compute_block_op_norm
    % maximum eigenvalue of operator
    for q = 1:R
        No = size(W{1}, 1);
        Tw_ = spalloc(size(T{q}, 1), No, size(T{q}, 2) * 16);
        Tw_(:, W{q}) = T{q};
        fprintf('\nComputing operator norm: block %i \n', q)
        evl_blocks(q) = op_norm(@(x) sqrt(aW{q}) .* (Tw_ * A(x)), @(x) At(Tw_' * (sqrt(aW{q}) .* x)), [Ny, Nx], 1e-6, 200, verbosity);
        clear Tw_;
    end
end

%% save data
if save_data_on_disk == 1
    fprintf('Saving new data ... \n');
    
    if save_data_on_disk
        file_start = save_dataset_number;
        while exist(sprintf('%s%s_input_data.mat', save_path, int2str(file_start)), 'file') || ...
                exist(sprintf('$s%s_input_data_config.mat', save_path, int2str(file_start)), 'file')
            file_start = file_start + 1;
        end
        if file_start ~= save_dataset_number;
            fprintf('WARNING: Saving new data in file %d instead of %d \n\n', file_start, save_dataset_number);
        end
        
        %         save(sprintf('%s%s_input_data', save_path, int2str(file_start)), '-v7.3', ... % 'G', 'Gw', 'W'
        %             'y0f', 'yf', 'y', 'y0', 'nWw'); % do not save G, it is large and is faster to compute than saving
        %         save(sprintf('%s%s_input_data_config', save_path, int2str(file_start)), 'N', 'Ny', ...
        %             'Nx', 'uw', 'vw', 'u', 'v', 'uvidx', 'im', 'sampling_pattern', 'param_sampling', 'input_snr', 'input_snr_v', 'image_file_name', 'num_tests', 'use_real_visibilities');
        
        if strcmp(sampling_pattern, 'file')
            for k = 1:num_tests
                vis = yf{k};
                %                 save(sprintf('%s%s_vis_data_t%s', save_path, int2str(file_start), int2str(k)), 'vis');
                clear vis;
                
                vis = yf{k} - y0f{k};
                %                 save(sprintf('%s%s_vis_data_noise_t%s', save_path, int2str(file_start), int2str(k)), 'vis');
                clear vis;
            end
        end
    end
end


if gen_data == 2
    fprintf('Using data from workspace ... \n\n');
end

if free_memory
    % free memory of the whose Gw or the split G are not needed
    try
        if ~run_admm_bpcon && ~run_sdmm_bpcon && ~run_fb_nnls && ~run_krylov_nnls
            clear Gw;
            clear Tw;
        end
        
        if ~run_pdfb_bpcon_par_sim_rescaled_precond_wave_par && ~run_pdfb_bpcon_par_sim && ~run_pdfb_bpcon_par_sim_rescaled && ~run_pdfb_bpcon_par_sim_rescaled_precond &&...
                ~run_pdfb_bpcon_par_sim_rand_rescaled && ~run_pdfb_bpcon_par_sim_block_rand_rescaled && ...
                ~run_pdfb_bpcon_dist && ~run_pdfb_bpcon_dist_rescaled && ~run_admm_bpconpar && ~run_sdmm_bpconpar && ...
                ~run_pdfb_bpcon_par_sim_rescaled_rec_async && ...
                ~run_pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p && ...
                ~run_admm_bpconpar_wavepar
            clear G;
            clear T;
            clear W;
        end
    end
end


%% WTF算法

addpath('Others');
addpath('Utilities');
addpath('Toolbox\DFrT_Tensor');
addpath('Solvers');
% evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);
xishu=op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);
% 使用极化SARA生成的加噪的可见度函数
y22=y{1,1}{1,1};
% Dirtymap=dirtyn./xishu;%脏图
Dirtymap=real(At(Gw'*y22)./evl);
[row,column]=size(im);

% 初始化
gamma_WTF = 0.2;  iter_max_WTF = 60 ;
gamma = gamma_WTF ;
iter_max = iter_max_WTF ;
tol_error = 1e-6;
im_diff_err = 1;
iter = 1;
imgMasked = Dirtymap;
nLvl = 4 ;
psi = TPCTFs(nLvl,6,row,column);% 3,4,6
im_rec_last = imgMasked;
im_rec = imgMasked;
im_rec_hat = im_rec;
t_last = 1;

tstart_a = tic;
while (im_diff_err>=tol_error)&&(iter<=iter_max)
    temp1=A(im_rec_hat);
    temp2 = temp1(W{1,1});
    temp3=T{1,1}*temp2;

    im_rec=im_rec_hat + gamma*(At(Tw' *(y22 - temp3)))./xishu;%数据保真度更新
    coefs = psi*im_rec;% coefs为小波系数α，psi为紧框架D
    thr_coefs =  thr_bishrink(coefs);% BS双变量收缩
    lamuda=norm(im_rec,1)/(row*column);
    im_rec = psi'*thr_coefs; % thr_coefs为取阈值后的α

    t_now = (1 + sqrt(1 + 4*t_last*t_last))/2;% t(k+1)
    im_rec_hat = im_rec + (t_last - 1)/t_now*(im_rec - im_rec_last);%更新后的图像
    if iter>1
        im_diff_err = norm(im_rec_last(:) - im_rec(:),2)/norm(im_rec_last(:),2);%误差大小
    end
    % 更新im_rec x(k)、t_now t(k)
    im_rec_last = im_rec;
    t_last = t_now;
    iter = iter +1;
    soln_SARA=real(im_rec);
    soln_SARA(soln_SARA< 0) =0;
    x1=A(soln_SARA);
    x2 = x1(W{1,1});
    x3=T{1,1}*x2;
    temp1 =norm(y22 -x3,2);
    %SNR
    error = im - im_rec;
    SNR = 20 * log10(norm(im(:))/norm(error(:)));
    fprintf('IterNum = %3d, SNR = %.4f dB\n',iter,SNR);
end
tend = toc(tstart_a);
%显示参数----------------------------------------------------
im_rec(im_rec<0) =0;%确保图像中没有负值。
im_rec =real(im_rec);%去掉虚部，仅保留实部。
ay = y{1,1};
ys = A(im_rec);
residual = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn_WTF = real(residual)/evl;
Dirtymap=real(At(Gw'*y22)./evl);
Sigma_WTF = norm(residualn_WTF(:)) ./ norm(Dirtymap(:));

fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('**************************************************************\n');
fprintf('**************************************************************\n');
fprintf('Final_Iter = %3d\n' , iter-1);
fprintf('Time = %.4f s\n', tend);

%计算指标-----------------------------------------------
%SNR
error = im - im_rec;
SNR = 20 * log10(norm(im(:))/norm(error(:)));
fprintf('SNR = %.4f dB\n', SNR);
fprintf('Sigma_ALSB: %.8f dB\n', Sigma_WTF);

%显示重构图----------------------------------------------------
figure; imagesc(real(log10((im_rec))));set(gca,'YDir','normal');axis off;
colorbar, axis image;
caxis([f_log_min f_log_max]);
title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

%显示误差图--------------------------------------------------------
error_WTF = im_rec - im ;
figure; imagesc(real(((error_WTF))));set(gca,'YDir','normal');axis off;
colorbar, axis image;
caxis([-f_err f_err]);
title('Error map');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

%显示残差图--------------------------------------------------------
figure; imagesc(real(((residualn_WTF))));set(gca,'YDir','normal');axis off;
colorbar, axis image;
caxis([-f_err f_err]);
title('residual map');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);


%% ALSB算法

addpath('Others');
addpath('Utilities');
addpath('Toolbox\DFrT_Tensor');
addpath('Solvers');
addpath('ksvdbox13');

% evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);
xishu = op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);
y22 = y{1,1}{1,1}; 
% Dirtymap=dirtyn./xishu;
Dirtymap = real(At(Gw'*y22)./evl);
[row,column] = size(im);

% 初始化
iter_max = 10;
tol_error = 1e-6;
im_diff_err = 1;
iter = 1;
imgMasked = Dirtymap; 
ALSB_Thr = 0.02;
mu = 0.002;
Step_length = 0.6;
Inloop = 20;
InitImg=imgMasked;

X = InitImg;
U = zeros(size(X)); 
B = zeros(size(X));

im_rec_last = imgMasked;
im_rec = imgMasked;
im_rec_hat = im_rec; 

tstart_a = tic;
while (im_diff_err>=tol_error)&&(iter<=iter_max)

    X_hat = X;
    R = X_hat-B;
    
    [X_bar,Dict] = ALSB_Solver(R,ALSB_Thr);    
    U = X_bar; % Dict*α

    for ii = 1:Inloop
        temp1 = A(im_rec); 
        temp2 = temp1(W{1,1});
        temp3 = T{1,1}*temp2;
        im_rec = (At(Tw' *(temp3 - y22)))./xishu;
        im_rec = real(im_rec);

        D = im_rec + mu*(X_hat - U - B);
        X = X_hat - D*Step_length;
        X(X<0) = 0;
        X = real(X);
        X_hat = X; 
        im_rec = X_hat;
    end
    B = B - (X - U);
    CurImg = X;
    
    %SNR
    error = im - CurImg;
    SNR = 20 * log10(norm(im(:))/norm(error(:)));
    fprintf('IterNum = %3d, SNR = %.4f dB\n',iter,SNR)    

    if iter>1
        im_diff_err = norm(im_rec_last(:) - im_rec(:),2)/norm(im_rec_last(:),2);%误差大小，norm计算向量的范数
    end
    im_rec_last = im_rec;
    iter = iter +1;   
end

RecImg = X;

tend = toc(tstart_a);
fprintf('Time: %.4f s\n', tend);

%数据处理----------------------------------------------------
RecImg(RecImg<0) =0;
RecImg =real(RecImg);

%计算SNR、Sigma-----------------------------------------------
error = RecImg - im ;%计算误差
SNR = 20 * log10(norm(im(:))/norm(real(error(:))));%计算SNR
fprintf('SNR: %.4f dB\n', SNR);

ay = y{1,1};
ys = A(RecImg);
residual = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn_ALSB = real(residual)/evl;
Dirtymap=real(At(Gw'*y22)./evl);
Sigma_ALSB = norm(residualn_ALSB(:)) ./ norm(Dirtymap(:));
fprintf('Sigma_ALSB: %.8f dB\n', Sigma_ALSB);

%显示重构图----------------------------------------------------
figure; imagesc(real(log10((RecImg))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([f_log_min f_log_max]);title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

%显示误差图--------------------------------------------------------
figure; imagesc(real(((error))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_err f_err]);title('Error map');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

%显示残差图--------------------------------------------------------
figure; imagesc(real(((residualn_ALSB))));set(gca,'YDir','normal');axis off;
colorbar, axis image; 
caxis([-f_res f_res]);
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

