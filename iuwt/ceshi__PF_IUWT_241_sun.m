%% 生成太阳圆盘原图
% % 目标场景
% clear;close all;clc;
% 
% % load('exampleimages.mat');
% % PSF=PSF(130:370,231:471);   % PSF最大值点250,351
% 
% Tc=241;
% jiaodu= 3 ;
% xi0 = linspace(-jiaodu,jiaodu,Tc);
% eta0 = linspace(-jiaodu,jiaodu,Tc);
% [xi,eta] = meshgrid(xi0,eta0);
% 
% TB_init= ones(size(xi));
% 
% % r=10;
% r=Tc*(0.5/(jiaodu*2));
% for x=1:Tc
%     for y=1:Tc
%        if (x-Tc/2)^2+(y-Tc/2)^2<=r^2                              
%            TB_init(x,y)=100000;
%        end
%     end
% end
% figure;pcolor(xi,eta,real(TB_init));shading interp;colorbar;title('TB init');
% Original=TB_init;
%%  用exampleimages内的脏束和脏图 测试PF和IUWT算法 
% clear;clc;close all;
% load exampleimages;
% 
% %Test for the PF algorithm
% fprintf('Starting PF based FISTA deconvolver');
% [PF_Model,PF_Residual]=FISTA_deconvolver_reweighted(Dirtymap,PSF,0.01,100,1);
% %Or Test for the IUWT algorithm
% fprintf('Starting IUWT based FISTA deconvolver');
% [IUWT_Model,IUWT_Residual]=FISTA_deconvolver_reweighted(Dirtymap,PSF,0.0001,100,1,1,4);
% 
% figure;pcolor(real(Original));shading interp;colorbar;title('Original');
% figure;pcolor(real(Dirtymap));shading interp;colorbar;title('Dirtymap');
% figure;pcolor(real(PF_Model));shading interp;colorbar;title('PF Model');
% figure;pcolor(real(IUWT_Model));shading interp;colorbar;title('IUWT Model');


%% 测试PF和IUWT   导入不同的脏束  原图为星云
% %导入matlab T1_psf生成的脏束
% % addpath('C:\Users\DELL\Desktop\L1\PSF');
% % load('B_241_3.mat');
% % PSF=B_sum;
% % PSF_max=max(max(PSF));
% % PSF=PSF./PSF_max;
% % PSF=real(PSF);
% % PSF(121,122)=0.9999;%有两个点值为1，设另一个为0.9999
% 
% %导入jupyter T2副本生成的脏束
% addpath('C:\Users\DELL\Desktop\yx\射电天文\网格化与加权\Imaging-Tutorial-master');
% A = importdata('T2-Bimage-3.csv');  %读取csv文件  脏束
% PSF=A;
% % PSF=zeros(size(A));PSF(118:126,118:126)=A(118:126,118:126); %取点扩散函数中间部分值
% 
% % 点扩散函数横截面的波形
% [len1,~]=size(PSF);
% x_hengjie=PSF((floor(len1/2)+1),:);
% [~,len2]=size(x_hengjie);
% xx=linspace(1,len2,len2);
% figure;plot(xx,real(x_hengjie));title('脏束横截面波形');
% 
% % 导入原图
% % addpath('C:\Users\DELL\Desktop\L1\原图');
% % TBimage= imread('7.jpg'); % 421*450
% % TBimage1=rgb2gray(TBimage);
% % TBimage1=double(TBimage1);
% % Tc=241;
% % Original=TBimage1(1:241,1:241);% 原图
% 
% % 与脏束卷积加噪声得脏图
% Dirtymap=conv2(Original,PSF);
% noise=0.01;%加噪声
% [mm,nn]=size(Dirtymap);
% Dirtymap2=Dirtymap(:);
% Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
% Dirtymap=reshape(Dirtymap2,mm,nn);
% Dirtymap=Dirtymap(121:361,121:361);% 脏图
% 
% %Test for the PF algorithm
% fprintf('Starting PF based FISTA deconvolver');
% [PF_Model,PF_Residual]=FISTA_deconvolver_reweighted(Dirtymap,PSF,0.01,100,1);
% %Or Test for the IUWT algorithm
% fprintf('Starting IUWT based FISTA deconvolver');
% [IUWT_Model,IUWT_Residual]=FISTA_deconvolver_reweighted(Dirtymap,PSF,0.0001,100,1,1,4);
% 
% figure;pcolor(xi,eta,real(Original));shading interp;colorbar;title('Original');
% figure;pcolor(xi,eta,real(PSF));shading interp;colorbar;title('Dirtybeam');
% figure;pcolor(xi,eta,real(Dirtymap));shading interp;colorbar;title('Dirtymap');
% figure;pcolor(xi,eta,real(PF_Model));shading interp;colorbar;title('PF Model');
% figure;pcolor(xi,eta,real(IUWT_Model));shading interp;colorbar;title('IUWT Model');

%% 测试PF和IUWT    导入不同的脏束  原图为太阳圆盘0.5度范围
% clear;clc;close all;
% %导入matlab T1_psf生成的脏束
% % addpath('C:\Users\DELL\Desktop\L1\PSF');
% % load('B_241_3.mat');
% % PSF=B_sum;
% % PSF_max=max(max(PSF));
% % PSF=PSF./PSF_max;
% % PSF=real(PSF);
% % PSF(121,122)=0.9999;%有两个点值为1，设另一个为0.9999
% 
% %导入jupyter T2副本生成的脏束
% addpath('C:\Users\DELL\Desktop\yx\射电天文\网格化与加权\Imaging-Tutorial-master');
% A = importdata('T2-Bimage-3.csv');  %读取csv文件  脏束
% PSF=A;
% % PSF=zeros(size(A));PSF(118:126,118:126)=A(118:126,118:126); %取点扩散函数中间部分值
% 
% % 点扩散函数横截面的波形
% [len1,~]=size(PSF);
% x_hengjie=PSF((floor(len1/2)+1),:);
% [~,len2]=size(x_hengjie);
% xx=linspace(1,len2,len2);
% figure;plot(xx,real(x_hengjie));title('脏束横截面波形');
% 
% % 原图  太阳圆盘0.5度范围
% addpath('C:\Users\DELL\Desktop\L1\原图');
% TBimage= imread('太阳.jpg'); % 610*610
% TBimage1=rgb2gray(TBimage);
% TBimage1=double(TBimage1);
% [m,n]=size(TBimage1);
% jiaodu= 5 ; % 圆弧角度
% k=tand(jiaodu);% 斜率 y=kx+b  b=n/2 
% temp=zeros(size(TBimage1));
% for ii=m/2:m
%     temp((n/2:floor((ii-n/2)*k+n/2)),ii)=1;
% end
% TBimage1=TBimage1.*temp;
% figure;pcolor(real(TBimage1));shading interp;colorbar;title('TBimage1');
% Original=TBimage1(185:425,305:545);% 原图  241*241
% figure;pcolor(real(Original));shading interp;colorbar;title('Original');
% 
% Dirtymap=conv2(Original,PSF);
% noise=0.01;%加噪声
% [mm,nn]=size(Dirtymap);
% Dirtymap2=Dirtymap(:);
% Dirtymap2 = Dirtymap2+sqrt(noise*max(Dirtymap2))*randn(size(Dirtymap2));
% Dirtymap=reshape(Dirtymap2,mm,nn);
% Dirtymap=Dirtymap(121:361,121:361);% 脏图
% 
% %Test for the PF algorithm
% fprintf('Starting PF based FISTA deconvolver');
% [PF_Model,PF_Residual]=FISTA_deconvolver_reweighted(Dirtymap,PSF,0.01,100,1);
% %Or Test for the IUWT algorithm
% fprintf('Starting IUWT based FISTA deconvolver');
% [IUWT_Model,IUWT_Residual]=FISTA_deconvolver_reweighted(Dirtymap,PSF,0.0001,100,1,1,4);
% 
% figure;pcolor(real(Original));shading interp;colorbar;title('Original');
% figure;pcolor(real(PSF));shading interp;colorbar;title('PSF');
% figure;pcolor(real(Dirtymap));shading interp;colorbar;title('Dirtymap');
% figure;pcolor(real(PF_Model));shading interp;colorbar;title('PF Model');
% figure;pcolor(real(IUWT_Model));shading interp;colorbar;title('IUWT Model');

%% 实际测量数据















