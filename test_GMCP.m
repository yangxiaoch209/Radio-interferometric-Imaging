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

importdata 3c353_gdth.fits ;
im=double(ans);
im=im/(max(max(im)));
[Nx,Ny]=size(im);
im=imresize(im,512./Nx);
% im(im<1e-4)=0;
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

%GMCP
[Ny, Nx]=size(im);
xishu = op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);
[row,column] = size(im);
hardt = @(z) max(real(z), 0);
soft = @(z, T) sign(z) .* max(abs(z)-T, 0);
y22=y{1}{1};
Dirtymap = zeros(row,column);
imgMasked = Dirtymap;
gamma = 0.5;
step = 1;
tol_error = 1e-6;
im_diff_err=1;
iter_max = 100;
iter = 1;
P = length(Psit);
im_rec = imgMasked;
im_v = imgMasked;
gmcp_alpha = cell(P, 1);
gmcp_v = cell(P, 1);
for k = 1:P
    gmcp_alpha{k} = zeros(size(Psit{k}(im_rec)));
    gmcp_v{k} = zeros(size(Psit{k}(im_v)));
end
util_create_pool(P);
tall=0;
tstart_gmcp = tic;
while (im_diff_err>=tol_error)&&(iter<=iter_max)
    %% L2 ball
    tall_time=tic;
    temp1 = gather(A(gpuArray(im_rec)));
    temp1 = temp1(W{1,1});
    temp1 = T{1,1}*temp1;
    temp1 = Tw'*(temp1-y22);
    temp1 = gather(At(gpuArray(temp1)))./xishu;
    
    temp2 = im_rec-im_v;
    temp2 = gather(A(gpuArray(temp2)));
    temp2 = temp2(W{1,1});
    temp2 = T{1,1}*temp2;
    temp2 = Tw'*temp2;
    temp2 = gather(At(gpuArray(temp2)))./xishu;
    %% soft_value
    lamuda = norm(im_rec,1)/(row*column);
    temp3 = lamuda*step;
    %% L1 ball
    for k = 1:P
        f(k) = parfeval(@run_par_waverec, 4, gmcp_alpha{k},gmcp_v{k},temp1,temp2,temp3, Psit{k}, Psi{k},step, gamma);
    end
    %% gather
    im_rec = zeros(size(im_rec));
    im_v = zeros(size(im_v));
    u1 = cell(P, 1);
    v1 = cell(P, 1);
    for k = 1:P
        [idx, im_rec_, im_v_,gmcp_alpha_,gmcp_v_] = fetchNext(f);
        u1{idx} = im_rec_;
        v1{idx} = im_v_;
        gmcp_alpha{idx} = gmcp_alpha_;
        gmcp_v{idx} = gmcp_v_;
        im_rec = im_rec + u1{idx};
        im_v = im_v + v1{idx};
    end
    im_rec = hardt(im_rec);
    im_v = hardt(im_v);
    
    if iter>1
        im_diff_err = norm(im_rec_last(:) - im_rec(:),2)/norm(im_rec_last(:),2);%误差大小
    end
    im_rec_last = im_rec;
    iter = iter +1;
end
tend = toc(tstart_gmcp);
if isempty(gcp('nocreate')) == 0
    delete(gcp);
end
gmcp_time = tend;

%---------------------------------------------
fprintf('tend = %.8f\n',tend);
%SNR
y22=y{1,1}{1,1};
im_rec = im_rec./max(max(im_rec));
im_rec =real(im_rec);
im_rec(im_rec<0)=0;
ay = {y22};
ys = A(im_rec);
residual = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn = real(residual)/evl;

error = im - im_rec;
SNR = 20 * log10(norm(im(:))/norm(error(:)));
fprintf('SNR = %.2f\n',SNR);

%sigma
Dirtymap=real(At(Gw'*y22)./evl);
sigma = norm(residualn(:)) ./ norm(Dirtymap(:));
fprintf('sigma = %.5f\n',sigma);

function [im_rec_, im_v_,gmcp_alpha_,gmcp_v_] = run_par_waverec(gmcp_alpha,gmcp_v,temp1,temp2,temp3, Psit, Psi,step,gamma)
temp1 = step*Psit(temp1);
temp2 = step*gamma*Psit(temp2);
gmcp_B = gmcp_alpha-temp1+temp2;
gmcp_C = gmcp_v+temp2;
gmcp_alpha = sign(gmcp_B) .* max(abs(gmcp_B)-temp3, 0);
gmcp_v = sign(gmcp_C) .* max(abs(gmcp_C)-temp3, 0);
im_rec_ = Psi(gmcp_alpha);
im_v_ = Psi(gmcp_v);
gmcp_alpha_ = gmcp_alpha;
gmcp_v_ = gmcp_v;
end
