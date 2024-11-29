%测试用的
clear;
clc;
close all;
addpath('C:\Users\DELL\Desktop\T_Array_Sim\function')
addpath('C:\Users\DELL\Desktop\L1\原图')
addpath('C:\Users\DELL\Desktop\L1\工作区')
addpath('C:\Users\DELL\Desktop\L1\function')
addpath 'C:\Users\DELL\Desktop\成像'
%% 圆形阵列
shape_radius=10; % 圆形半径
particle_x =[]; % 粒子x矩阵
particle_y =[];  % 粒子y矩阵
N= 157 ;  % 圆形均分数量
for k=1:N           
    x =shape_radius*cos(2*k*pi/N);
    y =shape_radius*sin(2*k*pi/N);
    particle_x =[particle_x x];
    particle_y =[particle_y y]; 
end
% plot(particle_x,particle_y,'r.');title('阵列布局');
% axis equal;
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

% figure;
% plot(u1,v1,'.r')
% xlabel('u');
% ylabel('v'); 
% title('u-v');
% axis equal;
% 目标场景
xi0 = linspace(-1,1,700);
eta0 = linspace(-1,1,700);
[xi,eta] = meshgrid(xi0,eta0);
% 天线方向图仿真
theta = real(asind(sqrt(xi.^2+eta.^2)));
Fn = cosd(theta).^3;
u=u1;v=v1;
Band=20*1e6;%带宽
Freq = 1.4e9;%中心频率


%% CLEAN
load('exampleimages.mat');
%%
% clear Original
% load('TB.mat');
% % figure;pcolor(TB);shading interp;colorbar; title('原图');
% TB=TB(1:500,151:850);
% Original=TB;
% Original=Original./25;
% figure;pcolor(Original);shading interp;colorbar; title('原图');
% clear Dirtymap
% Dirtymap=conv2(Original,PSF);%脏图
% noise=0.01;%加噪声
% [mm nn]=size(Dirtymap);
% Dirtymap2=Dirtymap(:);
% Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
% Dirtymap=reshape(Dirtymap2,mm,nn);
% Dirtymap=Dirtymap(250:749,351:1050);
% figure;pcolor(Dirtymap);shading interp;colorbar;title('脏图');
% caxis([0 1000]);
%%
% Dirtymap=Dirtymap/16; %测试
[m n]=find(PSF==max(max(PSF)));
center=[m n];
[m,n]=size(Dirtymap);
%Computng the UV mask with the psf
UV=fft2(circshift(PSF,1-center));
maxUV=max(max(abs(UV)));
indx=find(abs(UV/maxUV)>0.05);
invuv=zeros(m,n);
invuv(indx)=1./(UV(indx).^2);

%参数初始化
flag=1;
flag=gpuArray(flag);
gama=0.9; %循环因子，小于等于1
model=zeros(size(PSF));
model=gpuArray(model);
BD=PSF;%脏束
BD=gpuArray(BD);
Residual=zeros(m,n);%残差


%步骤一
% I_dirty=fft2(Dirtymap);%二维傅里叶变换，得初始脏束
temp_Dirtymap=Dirtymap;
Residual=zeros(m,n);%残差

times=0;
i=gpuArray(i);
Pmax=max(max(abs(temp_Dirtymap)));
Pmax=gpuArray(Pmax);

for i=1:500
    for j=1:700
        if BD(i,j)<0
            BD(i,j)=0;
        end
    end
end


while flag
    %步骤二
    % 脏图的绝对值峰值Pmax，出现位置tk
    Pmax=max(max((temp_Dirtymap)));
    [tx,ty]=find((temp_Dirtymap)==Pmax);
    [tx_BD,ty_BD]=find(BD==max(max((BD))));
%     if numel(tx)>1
%         ty=tx(2,1);
%         tx=tx(1,1);
%     end
   
    %步骤三
    %脏束放到寻找到的峰值点处的中心位置，更新脏图    B= circshift(A,K,m);K表示所移位数,m=1时表示列移位,m=2时表示行移位
    BD=circshift(BD,tx-tx_BD,1);%行移位
    BD=circshift(BD,ty-ty_BD,2);%列移位

    [txx,tyy]=find(BD==1);%调试用
    
    
    temp_Dirtymap=temp_Dirtymap-gama*BD*Pmax; %迭代更新    罗湾 乘以最大峰值Pmax，该点扩散函数Pmax＞＞1，迭代速度快，恢复效果比不加好
%     temp_Dirtymap=temp_Dirtymap-gama*BD; %迭代更新
    
    %步骤四
    %记录该峰值点源的幅度值以及位置，将其记录进空图像model中
%     model(tx,ty)=Pmax;
    model(tx,ty)=temp_Dirtymap(tx,ty);%记录减去脏束后该峰值点源的幅度值以及位置
    
    %步骤五
%     Residual=Dirtymap-temp_Dirtymap;
%     Residual=Dirtymap-real(ifft2(UV.*fft2(I_dirty)));%残图
%     wucha=sum(sum(temp_Dirtymap))
    times=times+1; %输出循环次数
    Pmax %输出最大值
    if    Pmax< 2 ||times>350000 %达到阈值
        flag=0;
    end
end

%% 对PSF脏束修正
%重构图像中CLEAN算法的研究
% b=PSF;

%方法一 运用窗函数修改“脏束”
b=PSF;
[tx,ty]=size(b);
win = hann(tx);
b=b.*win;


%方法二 取脏束的中心部分，减小波瓣带来的负面响应
% temp=zeros(500,700);
% temp(247:253,347:353)=0.025;
% b=b.*temp;


%方法三 直接拟合成高斯型的洁束
% A CLEAN type Deconvolution Algorithm该文，为了提高分辨率, 把迭代用的脏束由真实脏束变成了理想洁束
% X = -350 : 1 :349;
% Y = -249 : 1: 250;
% sigma=1;
% [X,Y]=meshgrid(X,Y);
% W=exp(-(X.^2+Y.^2)/sigma.^2);
% Z=W/2*pi*sigma.^2;
% Zmax=max(max(Z));
% Z=Z/Zmax;%最大值位置250,351
% mesh(X,Y,Z);
% title('权重分布图');
% 
% b=PSF; %最大值位置250,351
% b=b.*Z;
%% 结果
% 结果为model与PSF卷积＋剩余残图
I=conv2(model,b);
I=I(251:750,351:1050);
I=I+temp_Dirtymap;
% 结果为model＋剩余残图
% I=model+abs(temp_Dirtymap);


figure;pcolor(Original);shading interp;colorbar; title('原图');
% figure;pcolor(Dirtymap);shading interp;colorbar;title('脏图'); 
figure;pcolor(model);shading interp;colorbar;title('模型图'); 
figure;pcolor(I);shading interp;colorbar;title('结果图'); 

figure;pcolor(Original);shading interp;colorbar; caxis([0 10]);title('原图');
% figure;pcolor(Dirtymap);shading interp;colorbar; caxis([0 10]);title('脏图'); 
figure;pcolor(model);shading interp;colorbar; caxis([0 10]);title('模型图');
figure;pcolor(I);shading interp;colorbar; caxis([0 10]);title('结果图'); 


bias_Dirty=Dirtymap-Original;%bias误差
bias_Dirty=bias_Dirty(:);%拉成二维
RMSE_Dirty=sqrt(mean((real(bias_Dirty)).^2))

bias_model=model-Original;%bias误差
bias_model=bias_model(:);%拉成二维
RMSE_model=sqrt(mean((real(bias_model)).^2))

bias=I-Original;%bias误差
bias=bias(:);%拉成二维
RMSE=sqrt(mean((real(bias)).^2))


