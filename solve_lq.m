function [im_rec] = solve_lq(y,iter_max,gamma,Psitw,Psiw)
xishu=op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);%计算最大特征值


y22=y{1,1}{1,1};% 可见度函数
Dirtymap=At(Gw'*y22)./xishu;%脏图

[row,column]=size(im);

% 初始化
gamma = gamma  ; % 步长
iter_max  = iter_max ;%最大迭代次数
tol_error = 1e-6;%停止误差
im_diff_err=1;%初始为1
iter=1;
imgMasked=Dirtymap;
im_rec_last = imgMasked;
im_rec = imgMasked;
im_rec_hat = im_rec;
t_last = 1;


t_now = 1;
tor=0;
cnt = 0;
flag = 1;
L=1;
tstart_a = tic;
while (im_diff_err>=tol_error)&&(iter<=iter_max)
    
    im_rec_hat_old=im_rec_hat;
    
    temp1=A(im_rec_hat);
    temp2 = temp1(W{1,1});
    temp3=T{1,1}*temp2;
    im_rec=im_rec_hat + gamma*(At(Tw' *(y22 - temp3)))./xishu;% 更新图像
    im_rec(im_rec<0) =0;
    %%%%
    coef=Psitw(im_rec);% Psitw 稀疏变换 得到稀疏系数
    temp222=coef;
    lamda=norm(im_rec,1)/(row*column);% 正则化参数 
    coefs = shrinkage_Lq(coef, 0.8, lamda,L);  
    im_rec=Psiw(coefs);% Psiw 稀疏逆变换   得到图像
    
    t_last = t_now;
    t_now = (0.03+ sqrt(0.5 + r*t_last*t_last))/2;%更新tk
    a = min(1, (t_last-1) /t_now);
    
    im_rec_hat = im_rec + a *(im_rec - im_rec_last);% 更新后的图像
    
    if iter>1
        im_diff_err = norm(im_rec_last(:) - im_rec(:),2)/norm(im_rec_last(:),2);%误差大小
    end
    
    % 重启策略
    vk = (im_rec_hat_old(:)-im_rec(:))'*(im_rec(:)-im_rec_last(:));
    if vk >= tor
        cnt = cnt + 1;
        if cnt>=1
             t_now = (0.5 + sqrt(0.5 + r*t_last*t_last))/2;%更新tk
        else
            t_now = 1;% 重启
        end
        im_rec_hat = im_rec;% 重启
    end
 
    im_rec_last = im_rec;
    t_last = t_now;
    iter = iter +1;
end