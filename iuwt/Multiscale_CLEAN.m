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
clear Original
TBimage=importdata('太阳全日面图像.jpg');
TBimage1=rgb2gray(TBimage);
TBimage1=double(TBimage1);
figure;pcolor(TBimage1);shading interp;colorbar; % set(gca,'ydir','reverse'); %反转
Original=TBimage1;
Original=Original./25;
figure;pcolor(Original);shading interp;colorbar;
clear Dirtymap
Dirtymap=conv2(Original,PSF);%脏图
noise=0.01;%加噪声
[mm nn]=size(Dirtymap);
Dirtymap2=Dirtymap(:);
Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
Dirtymap=reshape(Dirtymap2,mm,nn);
Dirtymap=Dirtymap(250:749,351:1050);
figure;pcolor(Dirtymap);shading interp;colorbar;title('脏图');
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
%% 高塔
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
% load('Tower.mat');
% % figure;pcolor(Image);shading interp;colorbar; title('原图');
% Image=Image(1:500,151:850);
% Original=Image;
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
Residual=zeros(m,n);%残差
% I_dirty=fft2(Dirtymap);%二维傅里叶变换，得初始脏束
temp_Dirtymap=Dirtymap;
Residual=zeros(m,n);%残差

times=0;
% i=gpuArray(i);
Pmax=max(max(temp_Dirtymap));
% Pmax=gpuArray(Pmax);

% % 权重分布图 
% X = -350 : 1 :349;
% Y = -250 : 1: 249;
% sigma=1;
% [X,Y]=meshgrid(X,Y);
% W=exp(-(X.^2+Y.^2)/sigma.^2);%高斯
% Z=W/2*pi*sigma.^2;
% Zmax=max(max(Z));
% Z=Z/Zmax;
% mesh(X,Y,Z);
% title('权重分布图');

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

sigma_4 =8; % 方差大小 
ksigma_4 = ksigma(ksize, sigma_4, center,bias);
% figure; mesh(ksigma_4);
temp_Dirtymap_4=conv2(ksigma_4,temp_Dirtymap);
temp_Dirtymap_4=temp_Dirtymap_4(351:850,351:1050);
% figure;pcolor(temp_Dirtymap_4);shading interp;colorbar;
% max_4=max(max(temp_Dirtymap_4))

sigma_5 = 16; % 方差大小 
ksigma_5 = ksigma(ksize, sigma_5, center,bias);
% figure; mesh(ksigma_4);
temp_Dirtymap_5=conv2(ksigma_5,temp_Dirtymap);
temp_Dirtymap_5=temp_Dirtymap_5(351:850,351:1050);
% figure;pcolor(temp_Dirtymap_5);shading interp;colorbar;
% max_5=max(max(temp_Dirtymap_5))



% % 方法二：选择抛物线函数作为卷积核函数  该方法和脏图卷积后结果不对
% % 球面波函数
% lamda=6328e-10;                      %波长，单位：米
% k=2*pi/lamda;                        %波数
% x0=0.0001;                             %点光源的x坐标，单位：米
% y0=0.0001;                            %点光源的y坐标，单位：米
% z=0.3;                                %观察面到点光源的垂直距离，单位：米
% L=0.005;                              %观察面的尺寸，单位：米
% x=linspace(-L/2,L/2,700);y=x;        %构建x坐标和y坐标
% [x,y]=meshgrid(x,y);                 %构建二维坐标网格
% U1=exp(j*k*z).*exp(j*k.*((x-x0).^2+(y-y0).^2)/2/z);  %发散球面光波
% figure;pcolor(real(U1));shading interp;colorbar;
% Fai=real(U1);%取实数部分    
% % 抛物线函数
% % r=sqrt(x^2+y^2);
% % m=Fai.*(1-(r/a)^2);%抛物线函数
% % 若x,y为矩阵位置坐标
% a = 24 ;%a为当前尺度大小
% m=zeros(700,700);
% for i=1:700
%     for j=1:700
%         r=sqrt(i^2+j^2);
%         m(i,j)=Fai(i,j).*(1-(r/a)^2)/700;%抛物线函数
% %         if m(i,j)<0
% %             m(i,j)=0;
% %         end
%     end
% end
% figure;pcolor(real(m));shading interp;colorbar;
% m_max=max(max(real(m)))

%% 初始化 
% model=zeros(500,700);%空图像
% times=0;%循环次数
% BD=PSF;
% flag=1;
% gama= 0.5 ; %增益因子
% for i=1:500
%     for j=1:700
%         if BD(i,j)<0
%             BD(i,j)=0;
%         end
%     end
% end
% 
% % gpu加速
% temp_Dirtymap_1=gpuArray(temp_Dirtymap_1);
% temp_Dirtymap_2=gpuArray(temp_Dirtymap_2);
% temp_Dirtymap_3=gpuArray(temp_Dirtymap_3);
% temp_Dirtymap_4=gpuArray(temp_Dirtymap_4);
% temp_Dirtymap_5=gpuArray(temp_Dirtymap_5);
% flag=gpuArray(flag);
% model=gpuArray(model);
% BD=gpuArray(BD);
% times=gpuArray(times);
% tic %计时开始
% while flag
%     %% 步骤二 对于这些卷积脏图，找到全局最大峰值，并记下它的值、位置和尺度大小，储存在表里
%     max_temp_Dirtymap_1=max(max(temp_Dirtymap_1));%找到卷积脏图的峰值
%     max_temp_Dirtymap_2=max(max(temp_Dirtymap_2));
%     max_temp_Dirtymap_3=max(max(temp_Dirtymap_3));
%     max_temp_Dirtymap_4=max(max(temp_Dirtymap_4));
%     max_temp_Dirtymap_5=max(max(temp_Dirtymap_5));
%     temp=[max_temp_Dirtymap_1,max_temp_Dirtymap_2,max_temp_Dirtymap_3,max_temp_Dirtymap_4,max_temp_Dirtymap_5];
%     Pmax=max(temp);% 找到卷积脏图组中的全局最大峰值
%     if Pmax==max_temp_Dirtymap_1
%         max_temp_Dirtymap=temp_Dirtymap_1;
%     else
%         if Pmax==max_temp_Dirtymap_2
%             max_temp_Dirtymap=temp_Dirtymap_2;
%         else
%             if Pmax==max_temp_Dirtymap_3
%                 max_temp_Dirtymap=temp_Dirtymap_3;
%             else
%                 if Pmax==max_temp_Dirtymap_4
%                     max_temp_Dirtymap=temp_Dirtymap_4; % 卷积脏图组中的全局最大峰值对应的卷积脏图
%                 else
%                     if Pmax==max_temp_Dirtymap_5
%                         max_temp_Dirtymap=temp_Dirtymap_5;
%                     end
%                 end
%             end
%         end
%     end
%     [tx,ty]=find((max_temp_Dirtymap)==Pmax); %全局最大峰值的位置
%     model(tx,ty)=Pmax;%记录该峰值点源的幅度值以及位置，将其记录进空图像model中
% %     model(tx,ty)=model(tx,ty)+Pmax;
%     %% 步骤三  把脏束的中心移到该全局最大峰值的位置上，让脏束乘上一个增益因子，然后在所有的卷积脏图上减去这个脏束，得到每个尺度的残差图
%     [tx_BD,ty_BD]=find(BD==max(max((BD))));
%     BD=circshift(BD,real(tx-tx_BD),1);%行移位
%     BD=circshift(BD,real(ty-ty_BD),2);%列移位
% %     % 所有的卷积脏图上减去乘上增益因子的脏束
% % temp_Dirtymap_1=temp_Dirtymap_1-gama*BD*Pmax; %迭代更新
% % temp_Dirtymap_2=temp_Dirtymap_2-gama*BD*Pmax;
% % temp_Dirtymap_3=temp_Dirtymap_3-gama*BD*Pmax;
% % temp_Dirtymap_4=temp_Dirtymap_4-gama*BD*Pmax;
% % temp_Dirtymap_5=temp_Dirtymap_5-gama*BD*Pmax;
%     temp_Dirtymap_1=temp_Dirtymap_1-gama*BD; %迭代更新
%     temp_Dirtymap_2=temp_Dirtymap_2-gama*BD;
%     temp_Dirtymap_3=temp_Dirtymap_3-gama*BD;
%     temp_Dirtymap_4=temp_Dirtymap_4-gama*BD;
%     temp_Dirtymap_5=temp_Dirtymap_5-gama*BD;
%     times=times+1;  %循环次数
%     Pmax %输出最大值
%     if   Pmax < 1 || times> 350000  %达到阈值
%         flag=0;
%     end
% end
% toc % 计时
% %从gpu中取出
% model=gather(model);
% temp_Dirtymap_1=gather(temp_Dirtymap_1);
% temp_Dirtymap_2=gather(temp_Dirtymap_2);
% temp_Dirtymap_3=gather(temp_Dirtymap_3);
% temp_Dirtymap_4=gather(temp_Dirtymap_4);
% temp_Dirtymap_5=gather(temp_Dirtymap_5);
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
%         if temp_Dirtymap_4(i,j)<0
%             temp_Dirtymap_4(i,j)=0;
%         end
%         if temp_Dirtymap_5(i,j)<0
%             temp_Dirtymap_5(i,j)=0;
%         end
%     end
% end
%% 结果
% figure;pcolor(model);shading interp;colorbar; %模型图
% clean_1=conv2(model,ksigma_1);
% clean_1=clean_1(351:850,351:1050);
% clean_2=conv2(model,ksigma_2);
% clean_2=clean_2(351:850,351:1050);
% clean_3=conv2(model,ksigma_3);
% clean_3=clean_3(351:850,351:1050);
% clean_4=conv2(model,ksigma_4);
% clean_4=clean_4(351:850,351:1050);
% % 加剩余残图这一步是可选的  Parameterized CLEAN Deconvolution in Radio Synthesis Imaging
% % clean=clean_1+clean_2+clean_3+clean_4+temp_Dirtymap_1+temp_Dirtymap_2+temp_Dirtymap_3+temp_Dirtymap_4;%把所有的洁图和残图加起来
% clean=clean_1+clean_2+clean_3+clean_4;% 把所有的洁图加起来,最终的复原图
% figure;pcolor(Original);shading interp;colorbar;caxis([0 100]); title('原图');
% figure;pcolor(clean);shading interp;colorbar;title('复原图');
% 
% bias=clean-Original;%bias误差
% bias=bias(:);%拉成二维
% RMSE=sqrt(mean((real(bias)).^2))
%% 结果  model与不同尺度的核函数卷积＋所有的残图
% ksigma_11=ksigma_1(101:600,1:700);
% temp11=conv2(ksigma_11,PSF);
% temp11=temp11(251:750,351:1050);
% clean_1=conv2(model,temp11);
% clean_1=clean_1(251:750,351:1050);
% figure;pcolor(clean_1);shading interp;colorbar;
% 
% ksigma_22=ksigma_2(101:600,1:700);
% temp22=conv2(ksigma_22,PSF);
% temp22=temp22(251:750,351:1050);
% clean_2=conv2(model,temp22);
% clean_2=clean_2(251:750,351:1050);
% figure;pcolor(clean_2);shading interp;colorbar;
% 
% ksigma_33=ksigma_3(101:600,1:700);
% temp33=conv2(ksigma_33,PSF);
% temp33=temp33(251:750,351:1050);
% clean_3=conv2(model,temp33);
% clean_3=clean_3(251:750,351:1050);
% figure;pcolor(clean_3);shading interp;colorbar;
% 
% ksigma_44=ksigma_4(101:600,1:700);
% temp44=conv2(ksigma_44,PSF);
% temp44=temp44(251:750,351:1050);
% clean_4=conv2(model,temp44);
% clean_4=clean_4(251:750,351:1050);
% figure;pcolor(clean_4);shading interp;colorbar;
% 
% ksigma_55=ksigma_5(101:600,1:700);
% temp55=conv2(ksigma_55,PSF);
% temp55=temp55(251:750,351:1050);
% clean_5=conv2(model,temp55);
% clean_5=clean_5(251:750,351:1050);
% figure;pcolor(clean_5);shading interp;colorbar;
% 
% clean=clean_1+clean_2+clean_3+clean_4+clean_5;% 把所有的洁图加起来,最终的复原图
% clean2=clean_1+clean_2+clean_3+clean_4+clean_5+temp_Dirtymap_1+temp_Dirtymap_2+temp_Dirtymap_3+temp_Dirtymap_4+temp_Dirtymap_5;%把所有的洁图和残图加起来
% clean=clean./50;
% clean2=clean2./50;
% figure;pcolor(Original);shading interp;colorbar; %caxis([0 10]);title('原图');
% figure;pcolor(clean);shading interp;colorbar;%caxis([0 10]);title('复原图');
% figure;pcolor(clean2);shading interp;colorbar;%caxis([0 10]);title('复原图2');
% %复原图
% bias=clean-Original;%bias误差
% bias=bias(:);%拉成一维
% RMSE=sqrt(mean((real(bias)).^2))
% %复原图2
% bias2=clean2-Original;%bias误差
% bias2=bias2(:);%拉成一维
% RMSE2=sqrt(mean((real(bias2)).^2))
%%  多尺度clean    caj 基于压缩感知的太阳图像处理_罗湾
model=zeros(500,700);%空图像
times=0;%循环次数
BD=PSF;
flag=1;
gama= 1 ; %增益因子
for i=1:500
    for j=1:700
        if BD(i,j)<0
            BD(i,j)=0;
        end
    end
end

% gpu加速   gpu内存不足，比CPU慢
% temp_Dirtymap_1=gpuArray(temp_Dirtymap_1);
% temp_Dirtymap_2=gpuArray(temp_Dirtymap_2);
% temp_Dirtymap_3=gpuArray(temp_Dirtymap_3);
% temp_Dirtymap_4=gpuArray(temp_Dirtymap_4);
% temp_Dirtymap_5=gpuArray(temp_Dirtymap_5);
% ksigma_1=gpuArray(ksigma_1);
% ksigma_2=gpuArray(ksigma_2);
% ksigma_3=gpuArray(ksigma_3);
% ksigma_4=gpuArray(ksigma_4);
% ksigma_5=gpuArray(ksigma_5);
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
    max_temp_Dirtymap_4=max(max(temp_Dirtymap_4));
    max_temp_Dirtymap_5=max(max(temp_Dirtymap_5));
    temp=[max_temp_Dirtymap_1,max_temp_Dirtymap_2,max_temp_Dirtymap_3,max_temp_Dirtymap_4,max_temp_Dirtymap_5];
    Pmax=max(temp);% 找到卷积脏图组中的全局最大峰值
    if Pmax==max_temp_Dirtymap_1
        max_temp_Dirtymap=temp_Dirtymap_1;
        ksigma_temp=ksigma_1;
    else
        if Pmax==max_temp_Dirtymap_2
            max_temp_Dirtymap=temp_Dirtymap_2;
            ksigma_temp=ksigma_2;
        else
            if Pmax==max_temp_Dirtymap_3
                max_temp_Dirtymap=temp_Dirtymap_3;
                ksigma_temp=ksigma_3;
            else
                if Pmax==max_temp_Dirtymap_4
                    max_temp_Dirtymap=temp_Dirtymap_4; % 卷积脏图组中的全局最大峰值对应的卷积脏图
                    ksigma_temp=ksigma_4;
                else
                    if Pmax==max_temp_Dirtymap_5
                        max_temp_Dirtymap=temp_Dirtymap_5;
                        ksigma_temp=ksigma_5;
                    end
                end
            end
        end
    end
    [tx,ty]=find((max_temp_Dirtymap)==Pmax); %全局最大峰值的位置
%     model(tx,ty)=Pmax;%记录该峰值点源的幅度值以及位置，将其记录进空图像model中
    
    %% 步骤三  把脏束的中心移到该全局最大峰值的位置上，让脏束乘上一个增益因子，然后在所有的卷积脏图上减去这个脏束，得到每个尺度的残差图
    [tx_BD,ty_BD]=find(BD==max(max((BD))));
    BD=circshift(BD,real(tx-tx_BD),1);%行移位
    BD=circshift(BD,real(ty-ty_BD),2);%列移位
    % 尺度大小d和脏束BD进行卷积
    ksigma_temp=ksigma_temp(101:600,1:700);
    temp=conv2(ksigma_temp,BD);
    temp=temp(251:750,351:1050);
    model=model+gama.*temp;%更新model、
    
    % 所有的卷积脏图上减去乘上增益因子的脏束和尺度大小卷积
    temp_Dirtymap_1=temp_Dirtymap_1-gama.*temp.*Pmax; %迭代更新
    temp_Dirtymap_2=temp_Dirtymap_2-gama.*temp.*Pmax;
    temp_Dirtymap_3=temp_Dirtymap_3-gama.*temp.*Pmax;
    temp_Dirtymap_4=temp_Dirtymap_4-gama.*temp.*Pmax;
    temp_Dirtymap_5=temp_Dirtymap_5-gama.*temp.*Pmax;
%     temp_Dirtymap_1=temp_Dirtymap_1-gama*temp; %迭代更新
%     temp_Dirtymap_2=temp_Dirtymap_2-gama*temp;
%     temp_Dirtymap_3=temp_Dirtymap_3-gama*temp;
%     temp_Dirtymap_4=temp_Dirtymap_4-gama*temp;
%     temp_Dirtymap_5=temp_Dirtymap_5-gama*temp;
    times=times+1;  %循环次数
    Pmax %输出最大值
    if   Pmax < 3 || times> 350000  %达到阈值
        flag=0;
    end
end
toc % 计时
%从gpu中取出
% model=gather(model);
% temp_Dirtymap_1=gather(temp_Dirtymap_1);
% temp_Dirtymap_2=gather(temp_Dirtymap_2);
% temp_Dirtymap_3=gather(temp_Dirtymap_3);
% temp_Dirtymap_4=gather(temp_Dirtymap_4);
% temp_Dirtymap_5=gather(temp_Dirtymap_5);

% 取残图的正值
for i=1:500
    for j=1:700
        if temp_Dirtymap_1(i,j)<0
            temp_Dirtymap_1(i,j)=0;
        end
        if temp_Dirtymap_2(i,j)<0
            temp_Dirtymap_2(i,j)=0;
        end
        if temp_Dirtymap_3(i,j)<0
            temp_Dirtymap_3(i,j)=0;
        end
        if temp_Dirtymap_4(i,j)<0
            temp_Dirtymap_4(i,j)=0;
        end
        if temp_Dirtymap_5(i,j)<0
            temp_Dirtymap_5(i,j)=0;
        end
    end
end
%% 结果  存储结果model和洁束b卷积，并加上参与图，得到最终的洁图   caj基于压缩感知的太阳图像处理
%方法一 运用窗函数修改“脏束”
b=PSF;
[tx,ty]=size(b);
win = hann(tx);
b=b.*win;
% 结果为model与PSF卷积＋剩余残图
I=conv2(model,b);
I=I(251:750,351:1050);
% I=I+temp_Dirtymap_1+temp_Dirtymap_2+temp_Dirtymap_3+temp_Dirtymap_4+temp_Dirtymap_5;
% I=I./5;

% figure;pcolor(Original);shading interp;colorbar; caxis([0 10]);title('原图');
% figure;pcolor(Dirtymap);shading interp;colorbar;caxis([0 10]);title('脏图'); 
% figure;pcolor(model);shading interp;colorbar;caxis([0 10]);title('模型图'); 
figure;pcolor(I);shading interp;colorbar;caxis([0 10]);title('结果图'); 

bias=I-Original;%bias误差
bias=bias(:);%拉成二维
RMSE=sqrt(mean((real(bias)).^2))


