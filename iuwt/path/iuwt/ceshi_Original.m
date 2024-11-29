clear;
clc;
close all;
addpath('C:\Users\DELL\Desktop\T_Array_Sim\function')
addpath('C:\Users\DELL\Desktop\L1\原图')
addpath('C:\Users\DELL\Desktop\L1\工作区')
addpath('C:\Users\DELL\Desktop\L1\function')
%% 圆形阵列
shape_radius=10; % 圆形半径
particle_x =[]; % 粒子x矩阵
particle_y =[];  % 粒子y矩阵
N= 180 ;  % 圆形均分数量
for k=1:N           
    x =shape_radius*cos(2*k*pi/N);
    y =shape_radius*sin(2*k*pi/N);
    particle_x =[particle_x x];
    particle_y =[particle_y y]; 
end
plot(particle_x,particle_y,'r.');title('阵列布局');
axis equal;
CLAx=particle_x;
CLAy=particle_y;

%生成采样基线
[ak,aj] =pairs(length(CLAx));

%最短基线
du =0.6;
u0 = du*(CLAx(aj)-CLAx(ak));
v0 = du*(CLAy(aj)-CLAy(ak));

u1=[0;u0(:);-u0(:)];
v1=[0;v0(:);-v0(:)];

figure;
plot(u1,v1,'.r')
xlabel('u');
ylabel('v'); 
title('u-v');
axis equal;

% figure;
% plot(CLAx,CLAy,'r.');
% axis equal;
% xlabel('X');
% ylabel('Y');
% title('阵列布局');

%% 目标场景生成
load('Original.mat');
Original=Original(250:400,50:200);
[etaeta,xixi]=size(Original);
xi0 = linspace(-1,1,xixi);
eta0 = linspace(-1,1,etaeta);
[xi,eta] = meshgrid(xi0,eta0);
TB_init=Original;
dum = find(sqrt(xi.^2+eta.^2)<1);
xis = xi(dum);
etas = eta(dum);
TB = TB_init(dum);
TBs_init=nan(size(xi));
TBs_init(dum)=real(TB);
figure;pcolor(xi,eta,TBs_init);shading interp;colorbar; caxis([0 5]);%原始图像
%% 天线方向图仿真
theta = real(asind(sqrt(xi.^2+eta.^2)));
Fn = cosd(theta).^3;
Fns = Fn(dum);
u=u1;v=v1;
Band=20*1e6;%带宽
Freq = 1.4e9;%中心频率
%% G矩阵
% for ii=1:length(u) 
%   Gmatrix(ii,:)=Fns.*sinc(Band.*(u(ii)*xis+v(ii)*etas)./Freq).*exp(-1j*2*pi*(u(ii)*xis+v(ii)*etas))./sqrt(1-xis.^2-etas.^2);
% end
load('Gmatrix_180');
VF=Gmatrix*TB;
% temp=pinv(Gmatrix);
% TB2=temp*VF;
%% 加噪声
% noise=0.01;
% VF = VF+sqrt(noise*max(VF))*randn(size(VF));
%% 随机采样矩阵S m*M
% per=0.7;%数据使用量
% m=ceil(per*length(VF));
% S_eye=eye(length(VF));%生成单位阵
% temp1=(1:length(VF))';%位置
% temp2=temp1(randperm(length(temp1)));%随机顺序
% b=(temp2(1:(length(VF)-m)))';%要删除的行的标号
% S=S_eye;
% [row,col]=size(S_eye);   
% for i=1:length(b)
%     S(b(i),:)=ones(1,col)*NaN;            
% end
% S(isnan(S))=[];%删除完所有NaN后变成一列
% S=reshape(S,length(S)/col,col);%按照固定的列数重新排列,得到新矩阵
S=eye(length(VF));%百分百数据量
%% 离散余弦变换基D n*n
k = length(TB);          
D = sqrt(2 / k)*cos((0:k-1)'*((0:k-1)+0.5)*pi/k);
D(1,:) = D(1,:) / sqrt(2);%离散余弦变换基
% % T_dct = dct2(TB) ;    % matlab DCT函数DCT变换  idct逆变换
% % D=fft(eye(length(TB),length(TB)))/sqrt(length(TB));%离散傅里叶变换稀疏基
%% L1加权
f=S*VF;%测量数据
A=S*Gmatrix;
mu=10;%正则化系数
lambda=1;%拆分正则化系数
% Niter=10;%最大迭代次数
Niter=10;%最大迭代次数
a=1;%没有用到的参数
beta=300;%参数   60个阵元时23 180个阵元时300
%计算u = arg min ||u||_1+0.5*mu||Au-f||_2^2   u为TB
[u,mu,K,beta]=L1Adpt_tiduWWSplitBregmanIteration(f,A,mu,beta,Niter,a);
%%  计算T
% u2=D'*u;%逆变换
%无混叠视场
TB_resolve=nan(size(xi));
TB_resolve(dum)=real(u);
xi2=nan(size(xi(17:84,17:84)));
eta2=nan(size(eta(17:84,17:84)));
xi2(1:68,1:68)=xi(17:84,17:84);
eta2(1:68,1:68)=eta(17:84,17:84);
% figure;pcolor(xi2,eta2,TB_resolve(17:84,17:84));shading interp;colorbar; caxis([20 220]);
figure;pcolor(xi2,eta2,TB_resolve(17:84,17:84));shading interp;colorbar; caxis([0 5]);
bias=nan(size(xi));
bias=TB_resolve-TBs_init;%bias误差
bias_resolve=bias(17:84,17:84);
bias_resolve=bias_resolve(:);%拉成二维
RMSE_resolve=sqrt(mean((real(bias_resolve)).^2))
PSNR_resolve=20*(log10(255/RMSE_resolve));
% 无混叠原图
% figure;pcolor(xi2,eta2,TBs_init(17:84,17:84));shading interp;colorbar; caxis([20 220]);
%% --Tikhonov正则化 重构-------------
A = Gmatrix;b_bar=VF;
b = b_bar;
[U,s,V] = csvd (A);
figure
subplot (2,2,1);
lambda_l = l_curve (U,s,b); 
subplot (2,2,2);
k_l = l_curve (U,s,b,'tsvd'); 
subplot (2,2,3);
lambda_gcv = gcv (U,s,b); 
subplot (2,2,4);
k_gcv = gcv (U,s,b,'tsvd'); 
%---------Tikhonov正则化------------
T_tikh_l = tikhonov (U,s,V,b,lambda_l);
T_tikh_gcv = tikhonov (U,s,V,b,lambda_gcv);
% ------图像显示----------------
% TBs_tikh=nan(size(xi));
% TBs_tikh(dum)=real(T_tikh_gcv);
% figure;pcolor(xi,eta,TBs_tikh);shading interp; colorbar;caxis([20 220]);
% 求图像误差
% bias_tikh=TBs_tikh-TBs_init;
% 计算MAE误差
% MAE_tikh=mean(abs(real(T_tikh_gcv)-TB));
% 均方根误差
% RMSE_tikh=sqrt(1/length(dum)*sum((real(T_tikh_gcv)-TB).^2))

TB_tikh_gcv=nan(size(xi));
TB_tikh_gcv(dum)=real(T_tikh_gcv);
xi2=nan(size(xi(17:84,17:84)));
eta2=nan(size(eta(17:84,17:84)));
xi2(1:68,1:68)=xi(17:84,17:84);
eta2(1:68,1:68)=eta(17:84,17:84);
figure;pcolor(xi2,eta2,TB_tikh_gcv(17:84,17:84));shading interp;colorbar; caxis([20 220]);
bias=nan(size(xi));
bias=TB_tikh_gcv-TBs_init;%bias误差
bias_tikh_gcv=bias(17:84,17:84);
bias_tikh_gcv=bias_tikh_gcv(:);%拉成二维
RMSE_tikh_gcv=sqrt(mean((real(bias_tikh_gcv)).^2))
%% --------------TSVD正则化 重构-------------
%不加噪声反演不了
T_tsvd_l = tsvd (U,s,V,b,k_l);
T_tsvd_gcv = tsvd (U,s,V,b,k_gcv);
TBs_tsvd=nan(size(xi));
TBs_tsvd(dum)=real(T_tsvd_gcv);

TB_tsvd=nan(size(xi));
TB_tsvd(dum)=real(T_tsvd_gcv);
xi2=nan(size(xi(17:84,17:84)));
eta2=nan(size(eta(17:84,17:84)));
xi2(1:68,1:68)=xi(17:84,17:84);
eta2(1:68,1:68)=eta(17:84,17:84);
figure;pcolor(xi2,eta2,TB_tsvd(17:84,17:84));shading interp;colorbar; caxis([20 220]);
bias=nan(size(xi));
bias=TB_tsvd-TBs_init;%bias误差
bias_tsvd=bias(17:84,17:84);
bias_tsvd=bias_tsvd(:);%拉成二维
RMSE_tsvd=sqrt(mean((real(bias_tsvd)).^2))