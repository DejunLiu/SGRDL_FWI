%% Clear and close everything, open parallel pool
clear all
close all
clc
if isempty(gcp('nocreate')), parpool; end

%% Key parameters (designed to be edited quickly)
window_type = 'rectangle_outer';						% change to get rectangular window around receivers
fmin = 2;												% minimum frequency
fmax = 22;												% maximum frequency
nf = 12;												% number of frequencies
ctr = 1;												% contrast

%% Remaining parameters
% Ensures 10 grid points per wavelength
nx = ceil(384*10*fmax/122);									% number of gridpoints in x direction
ny = ceil(10*fmax);											% number of gridpoints in y direction
h=5e-3;
% number of gridpoints in y direction
dom = domain([0 (nx+1)*h 0 (ny+1)*h],[nx ny]);				% rectangular domain
wpml = 10*h;												% width of PML
freqs = linspace(fmin,fmax,nf);								% frequencies
ns = ceil(sqrt(nx*ny/(3*nf)));								% number of sources
nr = 3*ns;													% number of receivers
c_vec = [1 ctr+1];											% vector for colorbar
sigma =0;												    % noise level

maxit = 10;													% maximum number of LBFGS iterations
smoothness = 15;											% intensity of smoothing filter
source_info.type = 'hline';									% create a horizontal line of sources
source_info.bounds = [20*h (nx-20)*h];						% left and right endpoints of sources
source_info.height = 20*h;									% height of sources
sources = sources_and_receivers(ns,source_info);			% x and y locations of sources
receiver_info.type = 'hline';								% create a horizontal line of receivers
receiver_info.bounds = [20*h (nx-20)*h];								% left and right endpoints of receivers

receiver_info.height = 20*h;								% height of receivers
receivers = sources_and_receivers(nr,receiver_info);		% x and y locations of receivers
if strcmp(window_type,'all')								% decide whether or not to window
	window_info.type = 'all';								% take all points
else
	window_info.type = 'rectangle_outer';					% create a rectangular window
	window_info.bounds = [10*h (nx-10)*h 10*h (ny-10)*h];						% boundary of this window
end
[win_inds,W] = dom.window(window_info);						% indices of elements outside of window
pml_info.type = 'pml';										% pml info, for plotting purposes
pml_info.width = wpml;										% width of pml
[~,PML] = dom.window(pml_info);								% indicates whether a pixel is inside PML
v_true = marmousi(dom,ctr);					% true background velocity
c_true = dom.mat2vec(v_true);					% true background velocity
% initial model
v0 = @(zz,xx)v_true(1)+.75*max(zz-0.15,0);
v0 = v0(dom.Y,dom.X);
c0=dom.mat2vec(v0);
vmin=min(c_true);
vmax=max(c_true);

psnr = PSNR(v_true,v0);
rms  =   RMS(v_true,v0);
fprintf('Initial PSNR = %0.2f\n',psnr);
fprintf('Initial RMS = %0.2f\n',rms);





%pathes size
P_Size   	= 12;
S_Size      = 10;
%params.ori_imag=v_true;
params.T=16;
params.blocksize=P_Size;
params.stepsize=S_Size;
params.lambda=1e-3;
params.alpha=1e-3;


tic
for i=1:6
            
     if(i==1)
       c0 = dom.mat2vec(v0);	% initial estimate for background velocity
       u_l=0*c0;
       q_l=0*c0;
     else
       c0=dom.mat2vec(reconstructed_image);
       
     end

     [m,out] = adjoint_state_2d(dom,freqs,sources,receivers,window_info,c_true,c0,sigma,maxit,u_l,q_l);

     if sum(m<0)>0, warning('Negative values encountered in m!'), end
     index= find(m< 1./vmax^2);
     m(index) = 1./vmax^2;
     index= find(m> 1./vmin^2);
     m(index) = 1./vmin^2;
    [Data,D] = Buid_dic(dom.vec2mat(m+q_l), P_Size, S_Size);
     %% Deboising
     
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

     psnr =  PSNR(v_true,reconstructed_image);
     rms  =  RMS(v_true,reconstructed_image);

     fprintf('Initial PSNR = %0.2f\n',psnr);
     fprintf('Initial RMS = %0.2f\n',rms);
     
end


cpu_time = toc;
fprintf('Running time... %.2f (s).\n',cpu_time)

%% Post computation reconstrution and objective function
cpu_time = toc;
fprintf('Running time... %.2f (s).\n',cpu_time)


scrsz = get(0,'ScreenSize');


figure('Position',[1 scrsz(4)/4 scrsz(3)/2 scrsz(4)/4]);
imagesc(1e3*dom.x,1e3*dom.y,1e3*v_true);set(gca,'plotboxaspectratio',[7 2 2]);
caxis([1e3*vmin,1e3*vmax]);
xlabel('x [m]');ylabel('z [m]');
title('true model [m/s]');colorbar;


figure('Position',[1 scrsz(4)/4 scrsz(3)/2 scrsz(4)/4]);
imagesc(1e3*dom.x,1e3*dom.y,1e3*v0);set(gca,'plotboxaspectratio',[7 2 2]);
%caxis([1e3*vmin,1e3*vmax]);
xlabel('x [m]');ylabel('z [m]');
title('initial model [m/s]');colorbar;

% Reconstruction via FWI

% if sum(reconstructed_image<0)>0, warning('Negative values encountered in m!'), end
v_r=1e3*reconstructed_image;
figure('Position',[1 scrsz(4)/4 scrsz(3)/2 scrsz(4)/4]);
imagesc(1e3*dom.x,1e3*dom.y,v_r);set(gca,'plotboxaspectratio',[7 2 2]);
caxis([1e3*vmin,1e3*vmax]);
xlabel('x [m]');ylabel('z [m]');
title('inversion model [m/s]');colorbar;

% 
for offset=20:4:220
figure
plot(1e3*v0(:,offset),1e3*dom.y, 'k--', 1e3*v_true(:,offset),1e3*dom.y, 'b-.',v_r(:,offset), 1e3*dom.y, 'r');
set(gca,'plotboxaspectratio',[1 2 1]);
caxis([1e3*vmin,1e3*vmax]);
set(gcf,'unit','centimeters','position',[1 2 10 14]);
int1=floor(1e3*offset*h);
set(gca, 'YDir','reverse'); 
xlabel('v [m/s]');
ylabel('z [m]'); title(['x=',num2str(int1),'m']); legend('inital','true','SGRDL-FWI');
end


