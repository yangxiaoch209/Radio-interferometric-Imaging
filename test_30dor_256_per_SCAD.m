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
 importdata 30dor_256.fits
im=double(ans);
im=im./(max(max(im)));
[Nx,Ny]=size(im);
im=imresize(im,256./Nx);
[Nx,Ny]=size(im);N=Nx*Ny;

figure; imagesc(real(((im))));set(gca,'YDir','normal');axis off;
colorbar, axis image;%caxis([-2,0]);%title('Recovered image');0
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

%% 参数
input_snr =50; % noise level (on the measurements)



f_log_min=-2;
f_log_max=0;   
f_res=0.009;
f_err=0.06;

average_times=1;                                                                                                                                                                                                                                                                                                                                   

%%
per=0.5:0.1:0.5;
% per =0.4;
SNR_sum_SCAD=[];


TIME_sum_SCAD=[];

SNR_sum_ave_PD=zeros(1,length(per));

TIME_sum_ave_SCAD=zeros(1,length(per));
for iiiiiii=1:average_times
    
SNR_sum_SCAD=[];

TIME_sum_SCAD=[];

for iiiii=1:1:length(per)
param_sampling.p=per(1,iiiii);

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
% load ('test_ska_240s.uvw.mat');
% uw=uvw(:,1);vw=uvw(:,2);
% % MAX=max(max(real(uvw)));
% % uw=uw./MAX*pi;
% % vw=vw./MAX*pi;
figure;plot(uw,vw,'.black');
xlabel('\eta)','FontName','Times New Roman','FontSize',16);
ylabel('\eta','FontName','Times New Roman','FontSize',16);
xlim([-4,4]);
ylim([-4 4]);
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
Maxiter=150; soft_threshold_values=1e-4; forward_descent_step_size=0.49;%200
% run_pdfb_bpcon_par_sim_rescaled_precon
% PDFB parameter structure sent to the algorithm
param_pdfb_precond.im = im; % original image, used to compute the SNR
param_pdfb_precond.verbose = verbosity; % print log or not
param_pdfb_precond.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb_precond.nu2 = evl_precond; % bound on the norm of the operator A*G
% param_pdfb_precond.gamma = 1e-6; % convergence parameter L1 (soft th parameter)
param_pdfb_precond.gamma = soft_threshold_values ; % convergence parameter L1 (soft th parameter)
% param_pdfb_precond.tau = 0.49; % forward descent step size
param_pdfb_precond.tau = forward_descent_step_size ; % forward descent step size
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

% results
% script_gen_figures;


% save data
% script_save_result_data;



asol = result_st.sol{i};
residual1 = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn_PPD = real(residual1)/evl;

soln_PPD = asol;%/max(asol(:));
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




%% SCAD
xishu=op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);%计算最大特征值


y22=y{1,1}{1,1};% 可见度函数
Dirtymap=At(Gw'*y22)./xishu;%脏图
[row,column]=size(im);
p = 1 ;
q = 1 ;
r = 4 ;

% 初始化
gamma = 1  ; % 步长
iter_max  =80 ;%最大迭代次数
tol_error = 1e-6;%停止误差
im_diff_err=1;%初始为1
iter=1;
imgMasked=Dirtymap;


% 软阈值算子
soft = @(z, T) sign(z) .* max(abs(z)-T, 0);

%稀疏平均算子
nlevel=4;
dwtmode('per');
[C,S]=wavedec2(im,nlevel,'db8');
ncoef=length(C);
[C1,S1]=wavedec2(im,nlevel,'db1');
ncoef1=length(C1);
[C2,S2]=wavedec2(im,nlevel,'db2');
ncoef2=length(C2);
[C3,S3]=wavedec2(im,nlevel,'db3');
ncoef3=length(C3);
[C4,S4]=wavedec2(im,nlevel,'db4');
ncoef4=length(C4);
[C5,S5]=wavedec2(im,nlevel,'db5');
ncoef5=length(C5);
[C6,S6]=wavedec2(im,nlevel,'db6');
ncoef6=length(C6);
[C7,S7]=wavedec2(im,nlevel,'db7');
ncoef7=length(C7);

Psit22 = @(x) [wavedec2(x,nlevel,'db1')'; wavedec2(x,nlevel,'db2')';wavedec2(x,nlevel,'db3')';...
    wavedec2(x,nlevel,'db4')'; wavedec2(x,nlevel,'db5')'; wavedec2(x,nlevel,'db6')';...
    wavedec2(x,nlevel,'db7')';wavedec2(x,nlevel,'db8')']/sqrt(8);

Psi22 = @(x) (waverec2(x(1:ncoef1),S1,'db1')+waverec2(x(ncoef1+1:ncoef1+ncoef2),S2,'db2')+...
    waverec2(x(ncoef1+ncoef2+1:ncoef1+ncoef2+ncoef3),S3,'db3')+...
    waverec2(x(ncoef1+ncoef2+ncoef3+1:ncoef1+ncoef2+ncoef3+ncoef4),S4,'db4')+...
    waverec2(x(ncoef1+ncoef2+ncoef3+ncoef4+1:ncoef1+ncoef2+ncoef3+ncoef4+ncoef5),S5,'db5')+...
    waverec2(x(ncoef1+ncoef2+ncoef3+ncoef4+ncoef5+1:ncoef1+ncoef2+ncoef3+ncoef4+ncoef5+ncoef6),S6,'db6')+...
    waverec2(x(ncoef1+ncoef2+ncoef3+ncoef4+ncoef5+ncoef6+1:ncoef1+ncoef2+ncoef3+ncoef4+ncoef5+ncoef6+ncoef7),S7,'db7')+...
    waverec2(x(ncoef1+ncoef2+ncoef3+ncoef4+ncoef5+ncoef6+ncoef7+1:ncoef1+ncoef2+ncoef3+ncoef4+ncoef5+ncoef6+ncoef7+ncoef),S,'db8'))/sqrt(8);

% 初始化
im_rec_last = imgMasked;
im_rec = imgMasked;
im_rec_hat = im_rec;
t_last = 1;
t_now = 1;
tor=0;
cnt = 0;
flag = 1;
tstart_a = tic;
alphaa =0;
step_old=1;
step=step_old;

while (im_diff_err>=tol_error)&&(iter<=iter_max)
    gamma=step;
    im_rec_hat_old=im_rec_hat;
    temp1=A(im_rec_hat);
    temp2 = temp1(W{1,1});
    temp3=T{1,1}*temp2;
    im_rec=im_rec_hat +gamma*(At(Tw' *(y22 - temp3)))./xishu;% 更新图像
    temp11=A(im_rec);
    temp22 = temp11(W{1,1});
    temp33=T{1,1}*temp22;
    coef=Psitw(im_rec);% Psitw 稀疏变换 得到稀疏系数
 
     canshu=3e-3; %
     lamuda=log10(canshu*((norm((y22 - temp3),2))/(norm(coef,1)))+1) ;
     coefs=shrinkage_scad(coef,lamuda);
     im_rec=Psiw(coefs);% Psiw 稀疏逆变换   得到图像
     nei =At(Tw' *(temp33 - temp3))./xishu;
     no=norm(im_rec-im_rec_hat).^2;
    if sum((dot( nei ,im_rec-im_rec_hat)))>(0.95/step_old* no)
        step =0.90*norm(im_rec-im_rec_hat)^2./sum((dot(nei,im_rec-im_rec_hat)));
    else
        step=(1+1/(iter)^1.1)*step_old;
    end
    t_last = t_now;
    t_now = (p + sqrt(q + r*(step_old/step)*t_last*t_last))/2;%更新tk
    a = min(1, (t_last-1)/t_now);
    
    im_rec_hat = im_rec + a *(im_rec - im_rec_last);% 更新后的图像
   
  
    if iter>1
        im_diff_err = norm(im_rec_last(:) - im_rec(:),2)/norm(im_rec_last(:),2);%误差大小
    end
    vk = (im_rec_hat_old(:)-im_rec(:))'*(im_rec(:)-im_rec_last(:));
    if vk >= tor
        cnt = cnt + 1;
        if cnt>=1
             t_now = (1 + sqrt(1 + 4*t_last*t_last))/2;%更新tk
        else
            t_now = 1;% 重启
        end
        im_rec_hat = im_rec;% 重启
    end
    im_rec_last = im_rec;
    t_last = t_now;
    iter = iter +1;
    step_old=step;
end
tend = toc(tstart_a);
TIME_SCAD=tend;

im_rec =real(im_rec);%重构图像
im_rec(im_rec<0)=0;
soln_SCAD= im_rec;
% im_rec=im_rec/(max(max(im_rec)));
figure; imagesc(real(log10((soln_SCAD))));set(gca,'YDir','normal');axis off;
colorbar, axis image;caxis([ -2 0]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
ys = A(soln_SCAD);
residual = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn_SCAD = real(residual)./evl;
figure; imagesc(real(((residualn_SCAD))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_res f_res]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

error_SCAD =real( im-soln_SCAD) ;
figure; imagesc(real(((error_SCAD))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_err f_err]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
SNR_SCAD = 20 * log10(norm(im(:))/norm(error_SARA(:)))


%%
% SNR_sum_PD=[SNR_sum_PD SNR_PD];
SNR_sum_SCAD=[SNR_sum_SCAD SNR_SCAD];

% TIME_sum_PD=[TIME_sum_PD TIME_PD];
TIME_sum_SCAD=[TIME_sum_SCAD TIME_SCAD];
close all;

end

SNR_sum_ave_SCAD=SNR_sum_ave_SCAD +SNR_sum_SCAD ;

TIME_sum_ave_SCAD=TIME_sum_ave_SCAD +TIME_sum_SCAD ;


end

% SNR_sum_ave_PD=SNR_sum_ave_PD / average_times ;
SNR_sum_ave_SCAD=SNR_sum_ave_SCAD / average_times ;

TIME_sum_ave_SCAD=TIME_sum_ave_SCAD / average_times ;

close all;