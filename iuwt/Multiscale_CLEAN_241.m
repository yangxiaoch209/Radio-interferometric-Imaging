%% Multi-scale CLEAN  多尺度CLEAN              PSF取101*101，原图取101*101
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
% load('exampleimages.mat');
% PSF=PSF(130:370,231:471);   % PSF最大值点250,351

% addpath('C:\Users\DELL\Desktop\L1\PSF');
% load('B_241.mat');
% PSF=B_sum;

% %导入jupyter T2副本生成的脏束
addpath('C:\Users\DELL\Desktop\yx\射电天文\网格化与加权\Imaging-Tutorial-master');
A = importdata('T2-Bimage-3.csv');  %读取csv文件  脏束
PSF=A;

% 点扩散函数横截面的波形
[len1,~]=size(PSF);
x_hengjie=PSF((floor(len1/2)+1),:);
[~,len2]=size(x_hengjie);
xx=linspace(1,len2,len2);
figure;plot(xx,real(x_hengjie));title('脏束横截面波形');
%% 太阳圆盘   0.5度范围
addpath('C:\Users\DELL\Desktop\L1\原图');
TBimage= imread('太阳.jpg'); % 610*610
TBimage1=rgb2gray(TBimage);
TBimage1=double(TBimage1);
[m,n]=size(TBimage1);
jiaodu= 5 ; % 圆弧角度
k=tand(jiaodu);% 斜率 y=kx+b  b=n/2 
temp=zeros(size(TBimage1));
for ii=m/2:m
    temp((n/2:floor((ii-n/2)*k+n/2)),ii)=1;
end
TBimage1=TBimage1.*temp;
figure;pcolor(real(TBimage1));shading interp;colorbar;title('TBimage1');
Tc=241;
Original=TBimage1(185:425,305:545);% 原图  241*241
figure;pcolor(real(Original));shading interp;colorbar;title('Original');
% 脏图
Dirtymap=conv2(Original,PSF);
noise=0.01;%加噪声
[mm,nn]=size(Dirtymap);
Dirtymap2=Dirtymap(:);
Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
Dirtymap=reshape(Dirtymap2,mm,nn);
Dirtymap=Dirtymap(121:361,121:361);% 脏图
%% 星云    7.jpg   取241*241
% clear Original
% TBimage= imread('7.jpg'); % 421*450
% TBimage1=rgb2gray(TBimage);
% TBimage1=double(TBimage1);
% Tc=241;
% Original=TBimage1(1:241,1:241);
% % figure;pcolor(TBimage1);shading interp;colorbar; % set(gca,'ydir','reverse'); %反转
% figure;pcolor(Original);shading interp;colorbar;
% clear Dirtymap
% Dirtymap=conv2(Original,PSF);%脏图
% noise=0.01;%加噪声
% [mm,nn]=size(Dirtymap);
% Dirtymap2=Dirtymap(:);
% Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
% Dirtymap=reshape(Dirtymap2,mm,nn);
% Dirtymap=Dirtymap(121:361,121:361);
% % figure;pcolor(real(Dirtymap));shading interp;colorbar;title('脏图');
%% 星云    101.jpg
% clear Original
% TBimage=importdata('101.jpg');
% TBimage1=rgb2gray(TBimage);
% TBimage1=double(TBimage1);
% Tc=101;
% % figure;pcolor(TBimage1);shading interp;colorbar; % set(gca,'ydir','reverse'); %反转
% Original=TBimage1;
% Original=Original.*2;
% % figure;pcolor(Original);shading interp;colorbar;
% clear Dirtymap
% Dirtymap=conv2(Original,PSF);%脏图
% noise=0.01;%加噪声
% [mm,nn]=size(Dirtymap);
% Dirtymap2=Dirtymap(:);
% Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
% Dirtymap=reshape(Dirtymap2,mm,nn);
% Dirtymap=Dirtymap(51:151,51:151);
% % figure;pcolor(Dirtymap);shading interp;colorbar;title('脏图');
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
% clear Original
% load('Tower1.mat'); % 128*128
% % figure;pcolor(Image);shading interp;colorbar; title('原图');
% Tc=101;
% Image=T2(1:101,28:128);  % 101*101
% Original=Image;
% % Original=Original./10;
% figure;pcolor(Original);shading interp;colorbar; title('原图');
% clear Dirtymap
% Dirtymap=conv2(Original,PSF);%脏图
% noise=0.01;%加噪声
% [mm,nn]=size(Dirtymap);
% Dirtymap2=Dirtymap(:);
% Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
% Dirtymap=reshape(Dirtymap2,mm,nn);
% Dirtymap=Dirtymap(51:151,51:151);
% figure;pcolor(Dirtymap);shading interp;colorbar;title('脏图');
%% 飞机  
% clear Original
% TBimage=importdata('plane_128×128.jpg');
% TBimage1=rgb2gray(TBimage);
% TBimage1=double(TBimage1);
% Original=TBimage1(1:100,1:100);
% % Original=Original./10;
% figure;pcolor(Original);shading interp;colorbar; title('原图');
% clear Dirtymap
% Dirtymap=conv2(Original,PSF);%脏图
% noise=0.01;%加噪声
% [mm,nn]=size(Dirtymap);
% Dirtymap2=Dirtymap(:);
% Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
% Dirtymap=reshape(Dirtymap2,mm,nn);
% Dirtymap=Dirtymap(51:150,51:150);
% figure;pcolor(Dirtymap);shading interp;colorbar;title('脏图');
%% 初始原图
% clear Original
% Original=zeros(100,100);
% Original(40:60,40:60)=30;
% figure;pcolor(Original);shading interp;colorbar; title('原图');
% clear Dirtymap
% Dirtymap=conv2(Original,PSF);%脏图
% noise=0.01;%加噪声
% [mm nn]=size(Dirtymap);
% Dirtymap2=Dirtymap(:);
% Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
% Dirtymap=reshape(Dirtymap2,mm,nn);
% Dirtymap=Dirtymap(51:150,51:150);
% figure;pcolor(Dirtymap);shading interp;colorbar; title('脏图');
%% CLEAN
% [m,n]=find(PSF==max(max(PSF)));
% center=[m n];
% [m,n]=size(Dirtymap);
% %Computng the UV mask with the psf
% UV=fft2(circshift(PSF,1-center));
% maxUV=max(max(abs(UV)));
% indx=find(abs(UV/maxUV)>0.05);
% invuv=zeros(m,n);
% invuv(indx)=1./(UV(indx).^2);

%参数初始化
% I_dirty=fft2(Dirtymap);%二维傅里叶变换，得初始脏束
[m,n]=size(Dirtymap);
temp_Dirtymap=Dirtymap;
Residual=zeros(m,n);%残差
Pmax=max(max(temp_Dirtymap));

%%  步骤一 选取大小不同的尺度，将脏图与不同大小的尺度图像进行卷积，得到卷积图像组 
% 多尺度 CLEAN 函数的卷积核函数

% 方法一：选取高斯核函数  
ksize =Tc; % 高斯大小   矩阵大小
center  = round(ksize/2);  
bias  = ksize*10/10;  

sigma_1 =0.1; % 方差大小 
ksigma_1 = ksigma(ksize, sigma_1, center,bias);
% figure; mesh(ksigma_1);
temp_Dirtymap_1=conv2(ksigma_1,temp_Dirtymap);
temp_Dirtymap_1=temp_Dirtymap_1(121:361,121:361);
% figure;pcolor(temp_Dirtymap_1);shading interp;colorbar;
% max_1=max(max(temp_Dirtymap_1))

sigma_2 =0.2; % 方差大小 
ksigma_2 = ksigma(ksize, sigma_2, center,bias);
% figure; mesh(ksigma_2);
temp_Dirtymap_2=conv2(ksigma_2,temp_Dirtymap);
temp_Dirtymap_2=temp_Dirtymap_2(121:361,121:361);
% figure;pcolor(temp_Dirtymap_2);shading interp;colorbar;
% max_2=max(max(temp_Dirtymap_2))

sigma_3 =0.4; % 方差大小 
ksigma_3 = ksigma(ksize, sigma_3, center,bias);
% figure; mesh(ksigma_3);
temp_Dirtymap_3=conv2(ksigma_3,temp_Dirtymap);
temp_Dirtymap_3=temp_Dirtymap_3(121:361,121:361);
% figure;pcolor(temp_Dirtymap_3);shading interp;colorbar;
% max_3=max(max(temp_Dirtymap_3))

sigma_4 =0.8; % 方差大小 
ksigma_4 = ksigma(ksize, sigma_4, center,bias);
% figure; mesh(ksigma_3);
temp_Dirtymap_4=conv2(ksigma_4,temp_Dirtymap);
temp_Dirtymap_4=temp_Dirtymap_4(121:361,121:361);
% figure;pcolor(temp_Dirtymap_3);shading interp;colorbar;
% max_3=max(max(temp_Dirtymap_3))

sigma_5 =1.6; % 方差大小 
ksigma_5 = ksigma(ksize, sigma_5, center,bias);
% figure; mesh(ksigma_3);
temp_Dirtymap_5=conv2(ksigma_5,temp_Dirtymap);
temp_Dirtymap_5=temp_Dirtymap_5(121:361,121:361);
% figure;pcolor(temp_Dirtymap_3);shading interp;colorbar;
% max_3=max(max(temp_Dirtymap_3))
%% 结果  model与不同尺度的核函数卷积＋所有的残图

%%  多尺度clean    caj 基于压缩感知的太阳图像处理_罗湾
% model=zeros(500,700);%空图像
% times=0;%循环次数
% BD=PSF;
% flag=1;
% gama= 0.5 ; %增益因子
% flag1=0;flag2=0;flag3=0; %权重标志位初始化
% 
% % for i=1:500
% %     for j=1:700
% %         if BD(i,j)<0
% %             BD(i,j)=0;
% %         end
% %     end
% % end
% 
% % gpu加速   gpu内存不足，比CPU慢
% % temp_Dirtymap_1=gpuArray(temp_Dirtymap_1);
% % temp_Dirtymap_2=gpuArray(temp_Dirtymap_2);
% % temp_Dirtymap_3=gpuArray(temp_Dirtymap_3);
% % ksigma_1=gpuArray(ksigma_1);
% % ksigma_2=gpuArray(ksigma_2);
% % ksigma_3=gpuArray(ksigma_3);
% % flag=gpuArray(flag);
% % gama=gpuArray(gama);
% % model=gpuArray(model);
% % BD=gpuArray(BD);
% % times=gpuArray(times);
% 
% tic %计时开始
% while flag
%     %% 步骤二 对于这些卷积脏图，找到全局最大峰值，并记下它的值、位置和尺度大小，储存在表里
%     max_temp_Dirtymap_1=max(max(temp_Dirtymap_1));%找到卷积脏图的峰值
%     max_temp_Dirtymap_2=max(max(temp_Dirtymap_2));
%     max_temp_Dirtymap_3=max(max(temp_Dirtymap_3));
%     temp=[max_temp_Dirtymap_1,max_temp_Dirtymap_2,max_temp_Dirtymap_3];
%     Pmax=max(temp);% 找到卷积脏图组中的全局最大峰值
%     if Pmax==max_temp_Dirtymap_1
%         max_temp_Dirtymap=temp_Dirtymap_1;
%         ksigma_temp=ksigma_1;
%         a=sigma_1;
%         flag1=1;
%     else
%         if Pmax==max_temp_Dirtymap_2
%             max_temp_Dirtymap=temp_Dirtymap_2;
%             ksigma_temp=ksigma_2;
%             a=sigma_2; 
%             flag2=1;
%         else
%             if Pmax==max_temp_Dirtymap_3
%                 max_temp_Dirtymap=temp_Dirtymap_3;% 卷积脏图组中的全局最大峰值对应的卷积脏图
%                 ksigma_temp=ksigma_3;
%                 a=sigma_3;
%                 flag3=1;
%             end
%         end
%     end
%     [tx,ty]=find((max_temp_Dirtymap)==Pmax); %全局最大峰值的位置
% %     model(tx,ty)=Pmax;%记录该峰值点源的幅度值以及位置，将其记录进空图像model中
%     
%     %% 步骤三  把脏束的中心移到该全局最大峰值的位置上，让脏束乘上一个增益因子，然后在所有的卷积脏图上减去这个脏束，得到每个尺度的残差图
%     [tx_BD,ty_BD]=find(BD==max(max((BD))));
%     BD=circshift(BD,real(tx-tx_BD),1);%行移位
%     BD=circshift(BD,real(ty-ty_BD),2);%列移位
%     
%     % 尺度大小d和脏束BD进行卷积   进行卷积的BD每次不一样，无法提出来加速
%     ksigma_temp=ksigma_temp(101:600,1:700);
%     temp=conv2(ksigma_temp,BD);
%     temp=temp(251:750,351:1050);
%     
%     model=model+gama.*temp;%更新model、
%     
%     % 所有的卷积脏图上减去乘上增益因子的脏束和尺度大小卷积
%     temp_Dirtymap_1=temp_Dirtymap_1-gama.*temp.*Pmax; %迭代更新
%     temp_Dirtymap_2=temp_Dirtymap_2-gama.*temp.*Pmax;
%     temp_Dirtymap_3=temp_Dirtymap_3-gama.*temp.*Pmax;
% %     temp_Dirtymap_1=temp_Dirtymap_1-gama*temp; %迭代更新
% %     temp_Dirtymap_2=temp_Dirtymap_2-gama*temp;
% %     temp_Dirtymap_3=temp_Dirtymap_3-gama*temp;
% 
%     % 更新权重
%     sum=[sigma_1,sigma_2,sigma_3];
%     amax=max(sum);
%     S=1-0.03*a/amax;
% 
% %峰值残图乘以权重
%     if flag1
%         temp_Dirtymap_1=temp_Dirtymap_1.*S;
%     else
%         if flag2
%             temp_Dirtymap_2=temp_Dirtymap_2.*S;
%         else
%             if flag3
%                 temp_Dirtymap_3=temp_Dirtymap_3.*S;
%             end
%         end
%     end
%     flag1=0;flag2=0;flag3=0;%标志位置零
% 
%     %所有残图乘以权重 
% %     aerfa=0.001;
% %     S1=1-aerfa*sigma_1/amax;
% %     temp_Dirtymap_1=temp_Dirtymap_1.*S1;
% %     S2=1-aerfa*sigma_2/amax;
% %     temp_Dirtymap_2=temp_Dirtymap_2.*S2;
% %     S3=1-aerfa*sigma_3/amax;
% %     temp_Dirtymap_3=temp_Dirtymap_3.*S3;
% %     times=times+1;  %循环次数
%     Pmax %输出最大值
%     if   Pmax < 1e-2 || times> 350000  %达到阈值
%         flag=0;
%     end
% end
% toc % 计时
% 
% %从gpu中取出
% % temp_Dirtymap_1=gather(temp_Dirtymap_1);
% % temp_Dirtymap_2=gather(temp_Dirtymap_2);
% % temp_Dirtymap_3=gather(temp_Dirtymap_3);
% % ksigma_1=gather(ksigma_1);
% % ksigma_2=gather(ksigma_2);
% % ksigma_3=gather(ksigma_3);
% % flag=gather(flag);
% % gama=gather(gama);d
% % model=gather(model);
% % BD=gather(BD);
% % times=gather(times);
% 
% % 取残图的正值
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
% %% 结果  存储结果model和洁束b卷积，并加上参与图，得到最终的洁图   caj基于压缩感知的太阳图像处理
% %方法一 运用窗函数修改“脏束”
% b=PSF;
% [tx,ty]=size(b);
% win = hann(tx);
% b=b.*win;
% I=conv2(model,b);
% I=I(251:750,351:1050);
% I=I+temp_Dirtymap_1+temp_Dirtymap_2+temp_Dirtymap_3;
% 
% figure;pcolor(Original);shading interp;colorbar; title('原图');%caxis([0 10]);
% figure;pcolor(Dirtymap);shading interp;colorbar; title('脏图');%caxis([0 10]);
% figure;pcolor(model);shading interp;colorbar;title('模型图'); %caxis([0 10]);
% figure;pcolor(temp_Dirtymap_1);shading interp;colorbar;title('残图1'); %caxis([0 10]);
% figure;pcolor(temp_Dirtymap_2);shading interp;colorbar;title('残图2'); %caxis([0 10]);
% figure;pcolor(temp_Dirtymap_3);shading interp;colorbar;title('残图3'); %caxis([0 10]);
% figure;pcolor(I);shading interp;colorbar;title('结果图'); %caxis([0 10]);
% 
% bias=I-Original;%bias误差
% bias=bias(:);%拉成二维
% RMSE=sqrt(mean((real(bias)).^2))

%% 多尺度clean    一种基于多尺度带通滤波的洁化算法与 GPU 实现
% 初始化
model_1=zeros(Tc,Tc);model_2=zeros(Tc,Tc);model_3=zeros(Tc,Tc);model_4=zeros(Tc,Tc);model_5=zeros(Tc,Tc);%空图像
times_1=0;times_2=0;times_3=0;times_4=0;times_5=0;%循环次数
BD_1=PSF;BD_2=PSF;BD_3=PSF;BD_4=PSF;BD_5=PSF;
flag_1=1;flag_2=1;flag_3=1;flag_4=1;flag_5=1;
gama= 0.5 ; %增益因子
times_break1=0;times_break2=0;times_break3=0;times_break4=0;times_break5=0;
sum=[sigma_1,sigma_2,sigma_3,sigma_4,sigma_5];
amax=max(sum);
Pmax1_min=inf;Pmax2_min=inf;Pmax3_min=inf;Pmax4_min=inf;Pmax5_min=inf;%初始值

beta=0.1;
yuzhi1= 0.1*max(max(temp_Dirtymap_1)) ; %迭代截止阈值
yuzhi2= 0.1*max(max(temp_Dirtymap_2)) ; %迭代截止阈值
yuzhi3= 0.1*max(max(temp_Dirtymap_3)) ; %迭代截止阈值
yuzhi4= 0.1*max(max(temp_Dirtymap_4)) ; %迭代截止阈值
yuzhi5= 0.1*max(max(temp_Dirtymap_5)) ; %迭代截止阈值

while flag_1
    %步骤二
    % 脏图的绝对值峰值Pmax，出现位置tk
    Pmax=max(max((temp_Dirtymap_1)));
    Pmax_last=Pmax;
    [tx,ty]=find((temp_Dirtymap_1)==Pmax);
    [tx_BD,ty_BD]=find(BD_1==max(max((BD_1))));
   
    %步骤三
    %脏束放到寻找到的峰值点处的中心位置，更新脏图    B= circshift(A,K,m);K表示所移位数,m=1时表示列移位,m=2时表示行移位
    BD_1=circshift(BD_1,tx-tx_BD,1);%行移位
    BD_1=circshift(BD_1,ty-ty_BD,2);%列移位
    
        % 尺度大小d和脏束BD进行卷积   进行卷积的BD每次不一样，无法提出来加速
    ksigma_temp=ksigma_1;
    temp=conv2(ksigma_temp,BD_1);
    temp=temp(121:361,121:361);
    
    temp_Dirtymap_1=temp_Dirtymap_1-gama*temp*Pmax; %迭代更新    罗湾 乘以最大峰值Pmax，该点扩散函数Pmax＞＞1，迭代速度快，恢复效果比不加好
%     temp_Dirtymap=temp_Dirtymap-gama*BD; %迭代更新
    Pmax=max(max((temp_Dirtymap_1)));
    
    % 记录最小的Pmax
    temp_vecter=[Pmax,Pmax_last,Pmax1_min];
    Pmax1_min=min(temp_vecter);
    if Pmax-Pmax1_min>Pmax1_min*0.1
        times_break1=times_break1+1;  %跳出循环
    end
    
    
%     aerfa=0.6;
%     S1=1-aerfa*sigma_1/amax;
%     temp_Dirtymap_1=temp_Dirtymap_1.*S1;
    
    
    %步骤四
    %记录该峰值点源的幅度值以及位置，将其记录进空图像model中
    model_1=model_1+gama.*temp;%更新model、
%     model_1(tx,ty)=Pmax;%记录该峰值点源的幅度值以及位置，将其记录进空图像model中

    times_1=times_1+1;  %循环次数 
    Pmax; %更新后的最大值
    if    Pmax< yuzhi1 || times_1 > 350000 || times_break1>10  %达到阈值
        flag_1=0;
    end
end

while flag_2
    %步骤二
    % 脏图的绝对值峰值Pmax，出现位置tk
    Pmax=max(max((temp_Dirtymap_2)));
    Pmax_last=Pmax;
    [tx,ty]=find((temp_Dirtymap_2)==Pmax);
    [tx_BD,ty_BD]=find(BD_2==max(max((BD_2))));
   
    %步骤三
    %脏束放到寻找到的峰值点处的中心位置，更新脏图    B= circshift(A,K,m);K表示所移位数,m=1时表示列移位,m=2时表示行移位
    BD_2=circshift(BD_2,tx-tx_BD,1);%行移位
    BD_2=circshift(BD_2,ty-ty_BD,2);%列移位
    
            % 尺度大小d和脏束BD进行卷积   进行卷积的BD每次不一样，无法提出来加速
    ksigma_temp=ksigma_2;
    temp=conv2(ksigma_temp,BD_2);
    temp=temp(121:361,121:361);
    
    temp_Dirtymap_2=temp_Dirtymap_2-gama*temp*Pmax; %迭代更新    罗湾 乘以最大峰值Pmax，该点扩散函数Pmax＞＞1，迭代速度快，恢复效果比不加好
%     temp_Dirtymap=temp_Dirtymap-gama*BD; %迭代更新
    Pmax=max(max((temp_Dirtymap_2)));

        % 记录最小的Pmax
    temp_vecter=[Pmax,Pmax_last,Pmax2_min];
    Pmax2_min=min(temp_vecter);
    if Pmax-Pmax2_min>Pmax2_min*0.1
        times_break2=times_break2+1;  %跳出循环
    end
    
    
%     aerfa=0.1;
%     S2=1-aerfa*sigma_2/amax;
%     temp_Dirtymap_2=temp_Dirtymap_2.*S2;
    
    %步骤四
    %记录该峰值点源的幅度值以及位置，将其记录进空图像model中
    model_2=model_2+gama.*temp;%更新model、
%     model_2(tx,ty)=Pmax;%记录该峰值点源的幅度值以及位置，将其记录进空图像model中
    
    times_2=times_2+1; %输出循环次数
    Pmax; %输出最大值
    if    Pmax< yuzhi2 ||times_2>350000 || times_break2>10  %达到阈值
        flag_2=0;
    end
end

while flag_3
    %步骤二
    % 脏图的绝对值峰值Pmax，出现位置tk
    Pmax=max(max((temp_Dirtymap_3)));
    Pmax_last=Pmax;
    [tx,ty]=find((temp_Dirtymap_3)==Pmax);
    [tx_BD,ty_BD]=find(BD_3==max(max((BD_3))));
   
    %步骤三
    %脏束放到寻找到的峰值点处的中心位置，更新脏图    B= circshift(A,K,m);K表示所移位数,m=1时表示列移位,m=2时表示行移位
    BD_3=circshift(BD_3,tx-tx_BD,1);%行移位
    BD_3=circshift(BD_3,ty-ty_BD,2);%列移位
    
            % 尺度大小d和脏束BD进行卷积   进行卷积的BD每次不一样，无法提出来加速
    ksigma_temp=ksigma_3;
    temp=conv2(ksigma_temp,BD_3);
    temp=temp(121:361,121:361);
    
    temp_Dirtymap_3=temp_Dirtymap_3-gama*temp*Pmax; %迭代更新    罗湾 乘以最大峰值Pmax，该点扩散函数Pmax＞＞1，迭代速度快，恢复效果比不加好
%     temp_Dirtymap=temp_Dirtymap-gama*BD; %迭代更新
    Pmax=max(max((temp_Dirtymap_3)));

            % 记录最小的Pmax
    temp_vecter=[Pmax,Pmax_last,Pmax3_min];
    Pmax3_min=min(temp_vecter);
    if Pmax-Pmax3_min>Pmax3_min*0.1
        times_break3=times_break3+1;  %跳出循环
    end
    
    
%     aerfa=0.1;
%     S3=1-aerfa*sigma_3/amax;
%     temp_Dirtymap_3=temp_Dirtymap_3.*S3;
    
    %步骤四
    %记录该峰值点源的幅度值以及位置，将其记录进空图像model中
    model_3=model_3+gama.*temp;%更新model、
%     model_3(tx,ty)=Pmax;%记录该峰值点源的幅度值以及位置，将其记录进空图像model中
    
    times_3=times_3+1; %输出循环次数
    Pmax; %输出最大值
    if    Pmax< yuzhi3 ||times_3>350000  || times_break3>10 %达到阈值
        flag_3=0;
    end
end

while flag_4
    %步骤二
    % 脏图的绝对值峰值Pmax，出现位置tk
    Pmax=max(max((temp_Dirtymap_4)));
    Pmax_last=Pmax;
    [tx,ty]=find((temp_Dirtymap_4)==Pmax);
    [tx_BD,ty_BD]=find(BD_4==max(max((BD_4))));
   
    %步骤三
    %脏束放到寻找到的峰值点处的中心位置，更新脏图    B= circshift(A,K,m);K表示所移位数,m=1时表示列移位,m=2时表示行移位
    BD_4=circshift(BD_4,tx-tx_BD,1);%行移位
    BD_4=circshift(BD_4,ty-ty_BD,2);%列移位
    
            % 尺度大小d和脏束BD进行卷积   进行卷积的BD每次不一样，无法提出来加速
    ksigma_temp=ksigma_4;
    temp=conv2(ksigma_temp,BD_4);
    temp=temp(121:361,121:361);
    
    temp_Dirtymap_4=temp_Dirtymap_4-gama*temp*Pmax; %迭代更新    罗湾 乘以最大峰值Pmax，该点扩散函数Pmax＞＞1，迭代速度快，恢复效果比不加好
%     temp_Dirtymap=temp_Dirtymap-gama*BD; %迭代更新
    Pmax=max(max((temp_Dirtymap_4)));

            % 记录最小的Pmax
    temp_vecter=[Pmax,Pmax_last,Pmax4_min];
    Pmax4_min=min(temp_vecter);
    if Pmax-Pmax4_min>Pmax4_min*0.1
        times_break4=times_break4+1;  %跳出循环
    end
    
    
%     aerfa=0.1;
%     S3=1-aerfa*sigma_3/amax;
%     temp_Dirtymap_3=temp_Dirtymap_3.*S3;
    
    %步骤四
    %记录该峰值点源的幅度值以及位置，将其记录进空图像model中
    model_4=model_4+gama.*temp;%更新model、
%     model_3(tx,ty)=Pmax;%记录该峰值点源的幅度值以及位置，将其记录进空图像model中
    
    times_4=times_4+1; %输出循环次数
    Pmax; %输出最大值
    if    Pmax< yuzhi4 ||times_4>350000  || times_break4 >10 %达到阈值
        flag_4=0;
    end
end

while flag_5
    %步骤二
    % 脏图的绝对值峰值Pmax，出现位置tk
    Pmax=max(max((temp_Dirtymap_5)));
    Pmax_last=Pmax;
    [tx,ty]=find((temp_Dirtymap_5)==Pmax);
    [tx_BD,ty_BD]=find(BD_5==max(max((BD_5))));
   
    %步骤三
    %脏束放到寻找到的峰值点处的中心位置，更新脏图    B= circshift(A,K,m);K表示所移位数,m=1时表示列移位,m=2时表示行移位
    BD_5=circshift(BD_5,tx-tx_BD,1);%行移位
    BD_5=circshift(BD_5,ty-ty_BD,2);%列移位
    
            % 尺度大小d和脏束BD进行卷积   进行卷积的BD每次不一样，无法提出来加速
    ksigma_temp=ksigma_5;
    temp=conv2(ksigma_temp,BD_5);
    temp=temp(121:361,121:361);
    
    temp_Dirtymap_5=temp_Dirtymap_5-gama*temp*Pmax; %迭代更新    罗湾 乘以最大峰值Pmax，该点扩散函数Pmax＞＞1，迭代速度快，恢复效果比不加好
%     temp_Dirtymap=temp_Dirtymap-gama*BD; %迭代更新
    Pmax=max(max((temp_Dirtymap_5)));

            % 记录最小的Pmax
    temp_vecter=[Pmax,Pmax_last,Pmax5_min];
    Pmax5_min=min(temp_vecter);
    if Pmax-Pmax5_min>Pmax5_min*0.1
        times_break5=times_break5+1;  %跳出循环
    end
    
    
%     aerfa=0.1;
%     S3=1-aerfa*sigma_3/amax;
%     temp_Dirtymap_3=temp_Dirtymap_3.*S3;
    
    %步骤四
    %记录该峰值点源的幅度值以及位置，将其记录进空图像model中
    model_5=model_5+gama.*temp;%更新model、
%     model_3(tx,ty)=Pmax;%记录该峰值点源的幅度值以及位置，将其记录进空图像model中
    
    times_5=times_5+1; %输出循环次数
    Pmax; %输出最大值
    if    Pmax< yuzhi5 ||times_5>350000  || times_break5>10 %达到阈值
        flag_5=0;
    end
end

%方法一 运用窗函数修改“脏束”
b=PSF;
[tx,ty]=size(b);
win = hann(tx);
b=b.*win;
b=(win').*b;

I_1=conv2(model_1,b);
I_1=I_1(121:361,121:361);

I_2=conv2(model_2,b);
I_2=I_2(121:361,121:361);

I_3=conv2(model_3,b);
I_3=I_3(121:361,121:361);

I_4=conv2(model_4,b);
I_4=I_4(121:361,121:361);

I_5=conv2(model_5,b);
I_5=I_5(121:361,121:361);

%  模型和卷积核函数卷积
% I_1=conv2(model_1,ksigma_1);
% I_1=I_1(51:150,51:150);
% I_2=conv2(model_2,ksigma_2);
% I_2=I_2(51:150,51:150);
% I_3=conv2(model_3,ksigma_3);
% I_3=I_3(51:150,51:150);
% I_4=conv2(model_4,ksigma_4);
% I_4=I_4(51:150,51:150);
% I_5=conv2(model_5,ksigma_5);
% I_5=I_5(51:150,51:150);

% I=model_1+model_2+model_3+model_4+model_5;
I=I_1+I_2+I_3+I_4+I_5;

%亮温取正值
for i=1:Tc
    for j=1:Tc
        if I(i,j)<0
            I(i,j)=0;
        end
    end
end



figure;pcolor(Original);shading interp;colorbar; title('原图');caxis([0 200]);
figure;pcolor(PSF);shading interp;colorbar; title('PSF');
figure;pcolor(real(Dirtymap));shading interp;colorbar; title('脏图');%caxis([0 10]);
figure;pcolor(real(model_1));shading interp;colorbar;title('模型图1'); %caxis([0 10]);
figure;pcolor(real(model_2));shading interp;colorbar;title('模型图2'); %caxis([0 10]);
figure;pcolor(real(model_3));shading interp;colorbar;title('模型图3'); %caxis([0 10]);
figure;pcolor(real(model_4));shading interp;colorbar;title('模型图4'); %caxis([0 10]);
figure;pcolor(real(model_5));shading interp;colorbar;title('模型图5'); %caxis([0 10]);
% figure;pcolor(temp_Dirtymap_1);shading interp;colorbar;title('残图1'); %caxis([0 10]);
% figure;pcolor(temp_Dirtymap_2);shading interp;colorbar;title('残图2'); %caxis([0 10]);
% figure;pcolor(temp_Dirtymap_3);shading interp;colorbar;title('残图3'); %caxis([0 10]);
% figure;pcolor(temp_Dirtymap_4);shading interp;colorbar;title('残图4'); %caxis([0 10]);
% figure;pcolor(temp_Dirtymap_5);shading interp;colorbar;title('残图5'); %caxis([0 10]);
figure;pcolor(real(I_1));shading interp;colorbar;title('结果图1'); %caxis([0 10]);
figure;pcolor(real(I_2));shading interp;colorbar;title('结果图2'); %caxis([0 10]);
figure;pcolor(real(I_3));shading interp;colorbar;title('结果图3'); %caxis([0 10]);
figure;pcolor(real(I_4));shading interp;colorbar;title('结果图4'); %caxis([0 10]);
figure;pcolor(real(I_5));shading interp;colorbar;title('结果图5'); %caxis([0 10]);
figure;pcolor(real(I));shading interp;colorbar;title('多尺度clean'); %caxis([0 200]);

bias=I-Original;%bias误差
bias=bias(:);%拉成二维
RMSE=sqrt(mean((real(bias)).^2))


