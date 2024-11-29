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
importdata W28.fits ;
im=double(ans);
im=im/(max(max(im)));
[Nx,Ny]=size(im);
im=imresize(im,1024./Nx);
[Nx,Ny]=size(im);N=Nx*Ny;

figure; imagesc(real(log10((im))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-4 0]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

%% 参数
input_snr = 25 ; % noise level (on the measurements)

Maxiter= 200 ;soft_threshold_values=1e-4; forward_descent_step_size=0.49; 
lambda_IUWT= 1e-6 ;  gamma_IUWT = 0.5 ;  iter_max_IUWT  =  100 ;


f_log_min=-3;
f_log_max=0;
f_res=0.0005;
f_err=0.005;

param_sampling.p=0.1;
%% generate the sampling pattern
param_sampling.N = N; % number of pixels in the image
param_sampling.Nox = ox*Nx; % number of pixels in the image
param_sampling.Noy = oy*Ny; % number of pixels in the image

[uw, vw, ~] = util_gen_sampling_pattern(sampling_pattern, param_sampling);


if use_symmetric_fourier_sampling
    uw = [uw; -uw];
    vw = [vw; -vw];
end

%%  SKA
load ('test_ska_240s.uvw.mat');
uw=uvw(:,1);vw=uvw(:,2);
% MAX=max(max(real(uvw)));
% uw=uw./MAX*pi;
% vw=vw./MAX*pi;
figure;plot(uw,vw,'.r');

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







%% PPD
% run_pdfb_bpcon_par_sim_rescaled_precond

% PDFB parameter structure sent to the algorithm
param_pdfb_precond.im = im; % original image, used to compute the SNR
param_pdfb_precond.verbose = verbosity; % print log or not
param_pdfb_precond.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb_precond.nu2 = evl_precond; % bound on the norm of the operator A*G
param_pdfb_precond.gamma = 5e-6; % convergence parameter L1 (soft th parameter)
param_pdfb_precond.tau = 0.49; % forward descent step size
param_pdfb_precond.rel_obj = 0;1e-3; % stopping criterion
% param_pdfb_precond.max_iter = 20; % max number of iterations
param_pdfb_precond.max_iter = Maxiter ;
param_pdfb_precond.lambda0 = 1; % relaxation step for primal update
param_pdfb_precond.lambda1 = 1; % relaxation step for L1 dual update
param_pdfb_precond.lambda2 = 1; % relaxation step for L2 dual update
param_pdfb_precond.sol_steps = [inf]; % saves images at the given iterations

param_pdfb_precond.use_proj_elipse_fb = 1;
param_pdfb_precond.elipse_proj_max_iter = 200;
param_pdfb_precond.elipse_proj_min_iter = 1;
param_pdfb_precond.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid

param_pdfb_precond.use_reweight_steps = 0;
param_pdfb_precond.use_reweight_eps = 0;
param_pdfb_precond.reweight_steps = [500:100:2000 inf];
param_pdfb_precond.reweight_rel_obj = 1e-5; % criterion for performing reweighting
param_pdfb_precond.reweight_min_steps_rel_obj = 50;
param_pdfb_precond.reweight_alpha = 0.01; % omega^(0)
param_pdfb_precond.reweight_alpha_ff = 0.5; % exponential decay factor for omega^(0)


result_st = [];
result_st.sol = cell(num_tests, 1);
result_st.L1_v = cell(num_tests, 1);
result_st.L1_vp = cell(num_tests, 1);
result_st.L2_v = cell(num_tests, 1);
result_st.L2_vp = cell(num_tests, 1);
result_st.time = cell(num_tests, 1);
result_st.delta_v = cell(num_tests, 1);
result_st.sol_v = cell(num_tests, 1);
result_st.snr_v = cell(num_tests, 1);
result_st.sol_reweight_v = cell(num_tests, 1);
result_st.no_sub_itr_v = cell(num_tests, 1);
result_st.xcorr_v = cell(num_tests, 1);

result_st.snr = cell(num_tests, 1);
result_st.sparsity = cell(num_tests, 1);
result_st.no_itr = cell(num_tests, 1);


for i = 1:num_tests
    % wavelet mode is a global variable which does not get transfered
    % to the workes; we need to set it manually for each worker
    dwtmode('per');
    
    fprintf('Test run %i:\n', i);
    
    tstart_a = tic;
    fprintf(' Running pdfb_bpcon_par_sim_rescaled_precond\n');
    [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
        result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, result_st.no_sub_itr_v{i}, ~, ~, ...
        result_st.sol_reweight_v{i}] ...
        = pdfb_bpcon_par_sim_rescaled_precond(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, aW, W, Psi, Psit, Psiw, Psitw, param_pdfb_precond);
    tend = toc(tstart_a);
    fprintf(' pdfb_bpcon_par_sim_rescaled_precond runtime: %ds\n\n', ceil(tend));
    
    result_st.time{i} = tend;
    if ~use_real_visibilities
        error = im - result_st.sol{i};
        result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
    end
    result_st.no_itr{i} = length(result_st.L1_v{i});
    
    wcoef = [];
    for q = 1:length(Psit)
        wcoef = [wcoef; Psit{q}(result_st.sol{i})];
    end
    result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
end

% results_prefix = 'pdfb_bpcon_par_sim_rescaled_precond';
% param_structure_name = 'param_pdfb_precond';

snr = 0;
sp = 0;
soln = zeros(Ny, Nx);
residualn = zeros(Ny, Nx);
dirtyn = zeros(Ny, Nx);
asol = result_st.sol{i};
ay = y{i};
aL1_v = result_st.L1_v{i};
aL2_v = result_st.L2_v{i};
aL1_vp = result_st.L1_vp{i};
aL2_vp = result_st.L2_vp{i};
% ay0 = y0{i};
aepsilon = epsilon{i};
aepsilonT = epsilonT{i};
R = length(ay);
ys = A(asol);



asol = result_st.sol{i};
residual = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn_PPD = real(residual)/evl;

soln_PPD = asol;

figure; imagesc(real(log10((soln_PPD))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([f_log_min f_log_max]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);


figure; imagesc(real(((residualn_PPD))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_res f_res]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

error_PDD = im - soln_PPD;
figure; imagesc(real(((error_PDD))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_err f_err]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

SNR_PPD = 20 * log10(norm(im(:))/norm(error_PDD(:)));
TIME_PPD=tend;


%% IUWT算法

addpath iuwt/


y22=y{1,1}{1,1};
xishu=op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);
Dirtymap=dirtyn./xishu;%脏图


N= 1 ;%大循环 

positiveflg=1;
waveletflg=1;
level=4;

W_IUWT = @(x) IUWT(x,level);% W means IUWT  小波变换
WT_IUWT = @(x) IUWT(x,-level);% WT  means inverse IUWT

[m,n]=size(Dirtymap);

im_temp=Dirtymap;% 初始图像为脏图
X_temp=W_IUWT(zeros(m,n));% 小波变换域下的系数
weight=ones(size(X_temp));  %权重


lambda=lambda_IUWT;
gamma = gamma_IUWT ;% 正则化参数γ
iter_max  =  iter_max_IUWT ; 

tol_error = 1e-6;%停止误差
im_diff_err=1;%图像当前误差，初始为1
iter=1;%当前迭代次数
im_rec_hat=Dirtymap;
im_rec_last=Dirtymap;
t_new=1;
X_old=X_temp;

IUWT_SNR_huitu=[ ];
tstart_a = tic;
while (im_diff_err>=tol_error)&&(iter<=iter_max)  
    t_old=t_new;
    temp1=A(im_rec_hat);
    temp2 = temp1(W{1,1});
    temp3=T{1,1}*temp2;
    im_rec=im_rec_hat + gamma*(At(Tw' *(y22 - temp3)))./xishu;%Fu为U*F  
    coefs=W_IUWT(im_rec);
    % Soft thresholding
    Shrink=abs(coefs)-lambda*weight;
    X_temp=sign(coefs).*((Shrink>0).*Shrink); % X_temp为τ收缩算子   论文中的公式（23）
    im_rec = WT_IUWT(X_temp); % thr_coefs为取阈值后的α
     %Updating t and X
    t_new=(1+sqrt(1+4*t_old^2))/2; % t_new为t（k+1）   论文中的公式（24）
    im_rec_hat=im_rec+(t_old-1)/t_new*(im_rec-im_rec_last);  % X为β（k+1），更新IUWT域中脏图的小波系数 维度为500*3500     论文中的公式（25）
     if iter>1
            im_diff_err = norm(im_rec_last(:) - im_rec(:),2)/norm(im_rec_last(:),2);%误差大小
     end
    % 更新im_rec x(k)、t_now t(k)
    im_rec_last = im_rec;
    t_old = t_new;
    weight=1./(abs(W_IUWT(im_rec))+.1);  % 权重
    iter = iter +1;
    temp=real(im_rec);
end
tend = toc(tstart_a); 
u_IUWT=real(im_rec);

soln_IUWT = u_IUWT;
figure; imagesc(real(log10((soln_IUWT))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([f_log_min f_log_max]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

ys = A(soln_IUWT);
residual = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn_IUWT = real(residual)/evl;
figure; imagesc(real(((residualn_IUWT))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_res f_res]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
error_IUWT = im - u_IUWT;
figure; imagesc(real(((error_IUWT))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_err f_err]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
SNR_IUWT = 20 * log10(norm(im(:))/norm(error_IUWT(:)));
TIME_IUWT =tend;

%%
% LPG
%lq
iter_max =100;
gamma=1;
tstart = tic;
im_rec =solve_lq(y,iter_max,gamma,Psitw,Psiw);
TIME_lq =tend;

soln_lq=real(im_rec);%重构图像
ys = A(soln_lq);

figure; imagesc(real(log10((soln_lq))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-3 0]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

residua_lq = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn_lq = real(residual_lq)/evl;
figure; imagesc(real(((residualn_lq))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-0.005 0.005]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

error_lq = (im - soln_lq);
figure; imagesc(((real(error_lq))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-0.005 0.005]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
SNR_SARA = 20 * log10(norm((im(:)))/norm(real(error_SARA(:))));
