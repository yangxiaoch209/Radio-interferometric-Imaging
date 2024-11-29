%% Multi-scale CLEAN  多尺度CLEAN
clear;
clc;
close all;
addpath('C:\Users\DELL\Desktop\T_Array_Sim\function')
addpath('C:\Users\DELL\Desktop\L1\原图')
addpath('C:\Users\DELL\Desktop\L1\工作区')
addpath('C:\Users\DELL\Desktop\L1\function')
addpath 'C:\Users\DELL\Desktop\L1\首次应用压缩感知到射电天文'
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


%% 载入原图、脏图、脏束
load('exampleimages.mat');
%% 太阳全日面图像   来源：基于多尺度clean算法的太阳射电图像处理_汪子涵
% clear Original
% TBimage=importdata('太阳全日面图像.jpg');
% TBimage1=rgb2gray(TBimage);
% TBimage1=double(TBimage1);
% % figure;pcolor(TBimage1);shading interp;colorbar; % set(gca,'ydir','reverse'); %反转
% Original=TBimage1;
% Original=Original./25;
% % figure;pcolor(Original);shading interp;colorbar;
% clear Dirtymap
% Dirtymap=conv2(Original,PSF);%脏图
% noise=0.01;%加噪声
% [mm nn]=size(Dirtymap);
% Dirtymap2=Dirtymap(:);
% Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
% Dirtymap=reshape(Dirtymap2,mm,nn);
% Dirtymap=Dirtymap(250:749,351:1050);
% % figure;pcolor(Dirtymap);shading interp;colorbar;title('脏图');
%% 太阳
% clear Original
% TBimage=importdata('太阳.jpg');
% TBimage1=rgb2gray(TBimage);
% TBimage1=double(TBimage1);
% TBimage1=TBimage1(55:554,6:605);
% temp=zeros(500,50);
% TBimage1=[temp TBimage1 temp];
% % figure;pcolor(TBimage1);shading interp;colorbar; %set(gca,'ydir','reverse');
% Original=TBimage1;
% Original=Original./25;
% figure;pcolor(Original);shading interp;colorbar; 
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
%% 星球
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
% % caxis([0 1000]);
%% 高塔
clear Original
load('Tower.mat');
% figure;pcolor(Image);shading interp;colorbar; title('原图');
Image=Image(1:500,151:850);
Original=Image;
Original=Original./250;
figure;pcolor(Original);shading interp;colorbar; title('原图');
clear Dirtymap
Dirtymap=conv2(Original,PSF);%脏图
noise=0.01;%加噪声
[mm,nn]=size(Dirtymap);
Dirtymap2=Dirtymap(:);
Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
Dirtymap=reshape(Dirtymap2,mm,nn);
Dirtymap=Dirtymap(250:749,351:1050);
figure;pcolor(Dirtymap);shading interp;colorbar;title('脏图');
%% 初始原图
% clear Original
% Original=zeros(500,700);
% Original(200:300,300:400)=30;
% Original(400:450,400:450)=30;
% figure;pcolor(Original);shading interp;colorbar; title('原图');
% clear Dirtymap
% Dirtymap=conv2(Original,PSF);%脏图
% noise=0.01;%加噪声
% [mm nn]=size(Dirtymap);
% Dirtymap2=Dirtymap(:);
% Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
% Dirtymap=reshape(Dirtymap2,mm,nn);
% Dirtymap=Dirtymap(250:749,351:1050);
% figure;pcolor(Dirtymap);shading interp;colorbar;caxis([0 10]); title('脏图');
%% CLEAN
[m,n]=find(PSF==max(max(PSF)));
center=[m n];
[m,n]=size(Dirtymap);
%Computng the UV mask with the psf
UV=fft2(circshift(PSF,1-center));
maxUV=max(max(abs(UV)));
indx=find(abs(UV/maxUV)>0.05);
invuv=zeros(m,n);
invuv(indx)=1./(UV(indx).^2);

%参数初始化
% I_dirty=fft2(Dirtymap);%二维傅里叶变换，得初始脏束
temp_Dirtymap=Dirtymap;
Residual=zeros(m,n);%残差
Pmax=max(max(temp_Dirtymap));

%%  步骤一 选取大小不同的尺度，将脏图与不同大小的尺度图像进行卷积，得到卷积图像组 
% 多尺度 CLEAN 函数的卷积核函数

% 方法一：选取高斯核函数  
ksize =700; % 高斯大小   矩阵大小
center  = round(ksize/2);  
bias  = ksize*10/10;  

sigma_1 =1; % 方差大小 
ksigma_1 = ksigma(ksize, sigma_1, center,bias);
% figure; mesh(ksigma_1);
temp_Dirtymap_1=conv2(ksigma_1,temp_Dirtymap);
temp_Dirtymap_1=temp_Dirtymap_1(351:850,351:1050);
% figure;pcolor(temp_Dirtymap_1);shading interp;colorbar;
% max_1=max(max(temp_Dirtymap_1))

sigma_2 =2; % 方差大小 
ksigma_2 = ksigma(ksize, sigma_2, center,bias);
% figure; mesh(ksigma_2);
temp_Dirtymap_2=conv2(ksigma_2,temp_Dirtymap);
temp_Dirtymap_2=temp_Dirtymap_2(351:850,351:1050);
% figure;pcolor(temp_Dirtymap_2);shading interp;colorbar;
% max_2=max(max(temp_Dirtymap_2))

sigma_3 =4; % 方差大小 
ksigma_3 = ksigma(ksize, sigma_3, center,bias);
% figure; mesh(ksigma_3);
temp_Dirtymap_3=conv2(ksigma_3,temp_Dirtymap);
temp_Dirtymap_3=temp_Dirtymap_3(351:850,351:1050);
% figure;pcolor(temp_Dirtymap_3);shading interp;colorbar;
% max_3=max(max(temp_Dirtymap_3))

%% 结果  model与不同尺度的核函数卷积＋所有的残图

%%  多尺度clean    caj 基于压缩感知的太阳图像处理_罗湾
model=zeros(500,700);%空图像
times=0;%循环次数
BD=PSF;
flag=1;
gama= 0.5 ; %增益因子
S=1; % 初始权重为1
% for i=1:500
%     for j=1:700
%         if BD(i,j)<0
%             BD(i,j)=0;
%         end
%     end
% end

% gpu加速   gpu内存不足，比CPU慢
% temp_Dirtymap_1=gpuArray(temp_Dirtymap_1);
% temp_Dirtymap_2=gpuArray(temp_Dirtymap_2);
% temp_Dirtymap_3=gpuArray(temp_Dirtymap_3);
% ksigma_1=gpuArray(ksigma_1);
% ksigma_2=gpuArray(ksigma_2);
% ksigma_3=gpuArray(ksigma_3);
% flag=gpuArray(flag);
% gama=gpuArray(gama);
% model=gpuArray(model);
% BD=gpuArray(BD);
% times=gpuArray(times);

tic %计时开始
while flag
    %% 步骤二 对于这些卷积脏图，找到全局最大峰值，并记下它的值、位置和尺度大小，储存在表里
    max_temp_Dirtymap_1=max(max(temp_Dirtymap_1));%找到卷积脏图的峰值
    max_temp_Dirtymap_2=max(max(temp_Dirtymap_2));
    max_temp_Dirtymap_3=max(max(temp_Dirtymap_3));
    temp=[max_temp_Dirtymap_1,max_temp_Dirtymap_2,max_temp_Dirtymap_3];
    Pmax=max(temp);% 找到卷积脏图组中的全局最大峰值
    if Pmax==max_temp_Dirtymap_1
        max_temp_Dirtymap=temp_Dirtymap_1;
        ksigma_temp=ksigma_1;
        a=sigma_1;
%         temp_Dirtymap_1=temp_Dirtymap_1.*S;
    else
        if Pmax==max_temp_Dirtymap_2
            max_temp_Dirtymap=temp_Dirtymap_2;
            ksigma_temp=ksigma_2;
            a=sigma_2;
%             temp_Dirtymap_2=temp_Dirtymap_2.*S;
        else
            if Pmax==max_temp_Dirtymap_3
                max_temp_Dirtymap=temp_Dirtymap_3;% 卷积脏图组中的全局最大峰值对应的卷积脏图
                ksigma_temp=ksigma_3;
                a=sigma_3;
%                 temp_Dirtymap_3=temp_Dirtymap_3.*S;
            end
        end
    end
    [tx,ty]=find((max_temp_Dirtymap)==Pmax); %全局最大峰值的位置
%     model(tx,ty)=Pmax;%记录该峰值点源的幅度值以及位置，将其记录进空图像model中
    
    %% 步骤三  把脏束的中心移到该全局最大峰值的位置上，让脏束乘上一个增益因子，然后在所有的卷积脏图上减去这个脏束，得到每个尺度的残差图
    [tx_BD,ty_BD]=find(BD==max(max((BD))));
    BD=circshift(BD,real(tx-tx_BD),1);%行移位
    BD=circshift(BD,real(ty-ty_BD),2);%列移位
    % 尺度大小d和脏束BD进行卷积   进行卷积的BD每次不一样，无法提出来加速
    ksigma_temp=ksigma_temp(101:600,1:700);
    temp=conv2(ksigma_temp,BD);
    temp=temp(251:750,351:1050);
    
    model=model+gama.*temp;%更新model、
    
    % 所有的卷积脏图上减去乘上增益因子的脏束和尺度大小卷积
    temp_Dirtymap_1=temp_Dirtymap_1-gama.*temp.*Pmax; %迭代更新
    temp_Dirtymap_2=temp_Dirtymap_2-gama.*temp.*Pmax;
    temp_Dirtymap_3=temp_Dirtymap_3-gama.*temp.*Pmax;
%     temp_Dirtymap_1=temp_Dirtymap_1-gama*temp; %迭代更新
%     temp_Dirtymap_2=temp_Dirtymap_2-gama*temp;
%     temp_Dirtymap_3=temp_Dirtymap_3-gama*temp;

    % 更新权重
%     sum=[sigma_1,sigma_2,sigma_3];
%     amax=max(sum);
%     S=1-0.6*a/amax;

%峰值残图乘以权重   这里错了，temp_Dirtymap_1 已经更新，和max_temp_Dirtymap不相等
%     if max_temp_Dirtymap==temp_Dirtymap_1
%         temp_Dirtymap_1=temp_Dirtymap_1.*S;
%     else
%         if max_temp_Dirtymap==temp_Dirtymap_2
%             temp_Dirtymap_2=temp_Dirtymap_2.*S;
%         else
%             if max_temp_Dirtymap==temp_Dirtymap_3
%                 temp_Dirtymap_3=temp_Dirtymap_3.*S;
%             end
%         end
%     end

%所有残图乘以权重
%     temp_Dirtymap_1=temp_Dirtymap_1.*S;
%     temp_Dirtymap_2=temp_Dirtymap_2.*S;
%     temp_Dirtymap_3=temp_Dirtymap_3.*S;

    times=times+1;  %循环次数
    Pmax %输出最大值
    if   Pmax < 1 || times> 350000  %达到阈值
        flag=0;
    end
end
toc % 计时

%从gpu中取出
% temp_Dirtymap_1=gather(temp_Dirtymap_1);
% temp_Dirtymap_2=gather(temp_Dirtymap_2);
% temp_Dirtymap_3=gather(temp_Dirtymap_3);
% ksigma_1=gather(ksigma_1);
% ksigma_2=gather(ksigma_2);
% ksigma_3=gather(ksigma_3);
% flag=gather(flag);
% gama=gather(gama);d
% model=gather(model);
% BD=gather(BD);
% times=gather(times);

% 取残图的正值
% for i=1:500
%     for j=1:700
%         if temp_Dirtymap_1(i,j)<0
%             temp_Dirtymap_1(i,j)=0;
%         end
%         if temp_Dirtymap_2(i,j)<0
%             temp_Dirtymap_2(i,j)=0;
%         end
%         if temp_Dirtymap_3(i,j)<0
%             temp_Dirtymap_3(i,j)=0;
%         end
%     end
% end
%% 结果  存储结果model和洁束b卷积，并加上参与图，得到最终的洁图   caj基于压缩感知的太阳图像处理
%方法一 运用窗函数修改“脏束”
b=PSF;
[tx,ty]=size(b);
win = hann(tx);
b=b.*win;
I=conv2(model,b);
I=I(251:750,351:1050);
I=I+temp_Dirtymap_1+temp_Dirtymap_2+temp_Dirtymap_3;

figure;pcolor(Original);shading interp;colorbar; title('原图');%caxis([0 10]);
figure;pcolor(Dirtymap);shading interp;colorbar; title('脏图');%caxis([0 10]);
figure;pcolor(model);shading interp;colorbar;title('模型图'); %caxis([0 10]);
figure;pcolor(temp_Dirtymap_1);shading interp;colorbar;title('残图1'); %caxis([0 10]);
figure;pcolor(temp_Dirtymap_2);shading interp;colorbar;title('残图2'); %caxis([0 10]);
figure;pcolor(temp_Dirtymap_3);shading interp;colorbar;title('残图3'); %caxis([0 10]);
figure;pcolor(I);shading interp;colorbar;title('结果图'); %caxis([0 10]);

bias=I-Original;%bias误差
bias=bias(:);%拉成二维
RMSE=sqrt(mean((real(bias)).^2))

%% 结果  不同尺度的脏束拟合得到的洁束分别卷积表里的值，得一系列的洁图，所有的洁图和残图加起来，就可得到最终的复原图
% 超宽带多频点太阳射电图像处理算法研究_潘冰清













