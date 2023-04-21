clear all;close all;clc
addpath(genpath('.'));
%% 参数设置
S       	= 20;                            		
P_Size   	= 8;
S_Size      = 6;
load seismic.mat;
C_img       = data{2}(1:256,1:256); 
%C_img       = data{2};
[h, w] 	  	= size(C_img);  						
N_img 	    = C_img + S*randn(h, w); 	
N_img(N_img > 255) = 255; 
N_img(N_img < 0)   = 0; 
PSNR_I=20*log10(255/std2(C_img-N_img));
%% 建立树和学习字典
[Data,D] = Buid_dic(C_img, P_Size, S_Size);
%% 去噪
params.ori_imag=C_img;
params.D=D;
params.T=6;
params.L=laplacian(Data','nn',6);
params.Y=Data;
params.blocksize=P_Size;
params.stepsize=S_Size;
params.alpha=0.2;
params.mu=0.16;
params.lambda=1e-5;
PSNR_O=20*log10(255/std2(C_img-N_img))
[X,O_img]=GRSC_ADMM(params);
%YD=params.D*X;
%O_img= col2imstep(YD, size(C_img) ,[P_Size, P_Size], [S_Size, S_Size]);
%O_img(O_img > 255) = 255; 
%O_img(O_img < 0)   = 0; 




toc;


%% 画图
figure;
subplot(1,3,1);imagesc(C_img);%colormap(seismic(2));
subplot(1,3,2);imagesc(N_img);%colormap(seismic(2));
subplot(1,3,3);imagesc(O_img);%colormap(seismic(2));
PSNR_O=20*log10(255/std2(C_img-O_img))
figure;plot(C_img(:,100));
hold on
plot(O_img(:,100),'r--');
error=norm(O_img-C_img)/norm(C_img)