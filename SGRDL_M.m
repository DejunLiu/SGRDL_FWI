% Illustrate usage of adjoint_state_2d and related files to recover the
% Marmousi velocity model.
% NOTICE: This script can be configured to save data from different runs.


%% Clear and close everything, open parallel pool
clear all
close all
clc
if isempty(gcp('nocreate')), parpool; end

%% Key parameters (designed to be edited quickly)
window_type = 'rectangle_outer';						% change to get rectangular window around receivers
save_data = 1;											% determine whether we save info
data_dir = '~D:\matlab\wavetools-master\marm_data/';	% directory to save data in (if save_data=1)
fmin = 2.9;												% minimum frequency
fmax = 26;												% maximum frequency
nf = 12;												% number of frequencies
ctr = 1;												% contrast

%% Remaining parameters
% Ensures 10 grid points per wavelength
nx = ceil(384*8*fmax/122);									% number of gridpoints in x direction8
ny = ceil(8*fmax);											% number of gridpoints in y direction8
h=4.8*1e-3;
dom = domain([0 (nx+1)*h 0 (ny+1)*h],[nx ny]);				% rectangular domain
% h=dom.hx;
wpml = 10*h;											    % width of PML
freqs = linspace(fmin,fmax,nf);								% frequencies
ns = ceil(sqrt(nx*ny/(3*nf)));								% number of sources, in 1:3 ratio to receivers
nr = 3*ns;													% number of receivers
c_vec = [1 ctr+1];											% vector for colorbar
SNR =15;%1.9e-2;										% noise level 5%
% SNR = 1/eps;												% noise level 0%
maxit = 10;													% maximum number of LBFGS iterations
smoothness =10;											    % intensity of smoothing filter
source_info.type = 'hline';									% create a horizontal line of sources
source_info.bounds = [10*h (nx-10)*h];						% left and right endpoints of sources
%source_info.bounds = [.1 2.9];								% left and right endpoints of sources
source_info.height = 0.85;%0.85;0.82		0.8*0.75		% height of sources
sources = sources_and_receivers(ns,source_info);			% x and y locations of sources
receiver_info.type = 'hline';								% create a horizontal line of receivers
receiver_info.bounds = [10*h (nx-10)*h];					% left and right endpoints of receivers
%receiver_info.bounds = [.1 2.9];							% left and right endpoints of receivers
receiver_info.height = 0.8;%0.8;0.75		0.75*	0.7		% height of receivers
receivers = sources_and_receivers(nr,receiver_info);		% x and y locations of receivers
%window_type='all';

if strcmp(window_type,'all')								% decide whether or not to window
	window_info.type = 'all';								% take all points
else
	window_info.type = 'rectangle_outer';					% create a rectangular window
	window_info.bounds = [6*h (nx-6)*h 0.75 0.90];	        %boundary of this window
end
[win_inds,W] = dom.window(window_info);						% indices of elements outside of window
pml_info.type = 'pml';										% pml info, for plotting purposes
pml_info.width = wpml;										% width of pml
[~,PML] = dom.window(pml_info);								% indicates whether a pixel is inside PML
v_true = marmousi(dom,ctr);					                % true background velocity
c_true = dom.mat2vec(v_true);				               	% true background velocity
v0  =mean(v_true,2);
v0     = repmat(v0,[1 nx]); 

vmin=min(c_true);
vmax=max(c_true);

fsim =  cal_fsim(v_true,v0,0,0);
ssim =  cal_ssim(v_true,v0,0,0);
psnr = PSNR(v_true,v0);
rms  =   RMS(v_true,v0);
fprintf('Initial fSIM= %0.3f\n',fsim);
fprintf('Initial SSIM= %0.3f\n',ssim);
fprintf('Initial PSNR = %0.3f\n',psnr);
fprintf('Initial RMS = %0.3f\n',rms);
%% Saving data to disk
% Determine how many computational experiments have been run
if save_data
	exper_num = 1; %#ok<*UNRCH>
	while exist(sprintf('%s/setup%d.fig',data_dir,exper_num),'file')
		exper_num = exper_num+1;
	end
end

%% Preliminary plots of experimental setup
figure

% True background velocity, sources and receivers
subplot(2,4,1:2)
dom.imagesc(v_true);
hold on
dom.plot(receivers(1,:),receivers(2,:),'rx')
dom.plot(sources(1,:),sources(2,:),'go')
hold off
axis image
title('Problem setup')

% Initial background velocity
subplot(2,4,3:4)
dom.imagesc(v0)
axis image
title('Initial background velocity')


% Sources, receivers and PML
subplot(2,4,5:6)
imagesc(v_true+ctr*PML)
hold on
dom.plot(receivers(1,:),receivers(2,:),'rx')
dom.plot(sources(1,:),sources(2,:),'go')
hold off
title('With PML')
axis image
caxis([c_vec(1) c_vec(2)])

% Sources, receivers and window
subplot(2,4,7:8)
imagesc(v_true+ctr*flipud(W))
hold on
dom.plot(receivers(1,:),receivers(2,:),'rx')
dom.plot(sources(1,:),sources(2,:),'go')
hold off
title('With window')
axis image
caxis([c_vec(1) c_vec(2)])

% Format and save to disk
set(gcf,'Position',[11 340 1362 462])


%pathes size
P_Size   	= 12;
S_Size      = 10;
params.T=12;
params.blocksize=P_Size;
params.stepsize=S_Size;
params.lambda=1.2;                                                                                                                                                                                                                                                                                                                                     ;%
params.alpha=2;

tic
for i=1:4
            
     if(i==1)
       c0 = dom.mat2vec(v0);
       u_l=0*c0;
       q_l=0*c0;
     else
       c0=dom.mat2vec(reconstructed_image);  
     end

     [m,out] = adjoint_state_2d(dom,freqs,sources,receivers,window_info,c_true,c0,SNR,maxit,u_l,q_l,v_true);
      if sum(m<0)>0, warning('Negative values encountered in m!'), end
     index= find(m< 1./vmax^2);
     m(index) = 1./vmax^2;
     index= find(m> 1./vmin^2);
     m(index) = 1./vmin^2;
     reconstructed_image_sart = dom.vec2mat(sqrt(1./abs(m)));
     fprintf('Initial ssim = %0.2f\n',cal_ssim(reconstructed_image_sart,v_true,0,0));
     fprintf('Initial PSNR = %0.2f\n',csnr(reconstructed_image_sart,v_true,0,0));
     fprintf('Initial RMS = %0.2f\n',RMS(reconstructed_image_sart,v_true));
     fprintf('Initial FSIM = %0.2f\n',cal_fsim(reconstructed_image_sart,v_true,0,0));
    [Data,D] = Buid_dic(dom.vec2mat(m+q_l), P_Size, S_Size);
     %% ШЅды
     
     params.D=D;
     params.L=laplacian(Data','nn',12);
     params.Y=Data;
     params.ori_imag=dom.vec2mat(m+q_l);
     [X_gamm,reconstructed_m]=GRSC_ADMM(params);
     u_l=dom.mat2vec(reconstructed_m);
     q_l=q_l+m-u_l;
     
     reconstructed_image=sqrt(1./abs(reconstructed_m));
     index= find(reconstructed_image > vmax);
     reconstructed_image(index) = vmax;
     index= find(reconstructed_image< vmin);
     reconstructed_image(index) =  vmin;
     fsim =  cal_fsim(v_true,reconstructed_image,0,0);
     ssim =  cal_ssim(v_true,reconstructed_image,0,0);
     psnr =  PSNR(v_true,reconstructed_image);
     rms  =  RMS(v_true,reconstructed_image);
     fprintf('Initial fSIM= %0.3f\n',fsim);
     fprintf('Initial SSIM= %0.3f\n',ssim);
     fprintf('Initial PSNR = %0.3f\n',psnr);
     fprintf('Initial RMS = %0.3f\n',rms);
     
end


cpu_time = toc;
fprintf('Running time... %.2f (s).\n',cpu_time)

%% Post computation reconstrution and objective function
cpu_time = toc;
fprintf('Running time... %.2f (s).\n',cpu_time)

scrsz = get(0,'ScreenSize');


% if sum(reconstructed_image<0)>0, warning('Negative values encountered in m!'), end
v_r=1e3*reconstructed_image;
figure('Position',[1 scrsz(4)/4 scrsz(3)/2 scrsz(4)/4]);
imagesc(1e3*dom.x,1e3*dom.y,v_r);set(gca,'plotboxaspectratio',[7 2 2]);
caxis([1e3*vmin,1e3*vmax]);
xlabel('x [m]');ylabel('z [m]');
title('SGRDL Reconstruction [m/s]');colorbar;




