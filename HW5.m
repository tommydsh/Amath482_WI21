%% Clean workspace
clear all; close all; clc
format shortG

%%
MCL = VideoReader('monte_carlo_low.mp4');
mcl_frames = read(MCL);
SDL = VideoReader('ski_drop_low.mp4');
sdl_frames = read(SDL);
numFrames_mcl = size(mcl_frames,4);
numFrames_sdl = size(sdl_frames,4);
t_mcl = linspace(0,MCL.Duration,numFrames_mcl); dt_mcl = t_mcl(2) - t_mcl(1);
t_sdl = linspace(0,SDL.Duration,numFrames_sdl); dt_sdl = t_sdl(2) - t_sdl(1);

%%
M_mcl = zeros(size(mcl_frames,1)*size(mcl_frames,2),numFrames_mcl);
M_sdl = zeros(size(sdl_frames,1)*size(sdl_frames,2),numFrames_sdl);
for i = 1:numFrames_mcl
   M_mcl(:,i) = reshape(im2double(rgb2gray(mcl_frames(:,:,:,i))),[size(mcl_frames,1)*size(mcl_frames,2),1]);
end
for i = 1:numFrames_sdl
   M_sdl(:,i) = reshape(im2double(rgb2gray(sdl_frames(:,:,:,i))),[size(sdl_frames,1)*size(sdl_frames,2),1]);
end

%% Monte-Carlo
M_mcl_1 = M_mcl(:,1:end-1);
M_mcl_2 = M_mcl(:,2:end);
[U_mcl,Sigma_mcl,V_mcl] = svd(M_mcl_1,'econ');
S_mcl = U_mcl'*M_mcl_2*V_mcl*diag(1./diag(Sigma_mcl));
[eV_mcl,D_mcl] = eig(S_mcl);
mu_mcl = diag(D_mcl);
omega_mcl = log(mu_mcl)/dt_mcl;
Phi_mcl = U_mcl*eV_mcl;
y0_mcl = Phi_mcl\M_mcl_1(:,1);

mode_mcl = 0;
mode_mcl_i = zeros(length(omega_mcl),1);
for i = 1:length(omega_mcl)
    if abs(real(omega_mcl(i))) < 0.01
        mode_mcl = mode_mcl+1;
        mode_mcl_i(mode_mcl) = i;
    end
end
mode_mcl_i = mode_mcl_i(1:mode_mcl);

y0_mcl_low = y0_mcl(mode_mcl_i);
omega_mcl_low = omega_mcl(mode_mcl_i);
Phi_mcl_low = Phi_mcl(:,mode_mcl_i);

u_modes_mcl_low = zeros(length(y0_mcl_low),numFrames_mcl);
for iter = 1:numFrames_mcl
   u_modes_mcl_low(:,iter) = y0_mcl_low.*exp(omega_mcl_low*t_mcl(iter)); 
end
u_dmd_mcl_low = Phi_mcl_low*u_modes_mcl_low;

u_dmd_mcl_sparse = M_mcl - abs(u_dmd_mcl_low);
u_mcl_R = zeros(size(u_dmd_mcl_sparse,1),size(u_dmd_mcl_sparse,2));
for i = 1:size(u_dmd_mcl_sparse,1)
    for j = 1:size(u_dmd_mcl_sparse,2)
        if u_dmd_mcl_sparse(i,j) < 0
            u_mcl_R(i,j) = u_dmd_mcl_sparse(i,j);
        end
    end
end
u_dmd_mcl_low_final = u_mcl_R+abs(u_dmd_mcl_low);
u_dmd_mcl_sparse_final = u_dmd_mcl_sparse-u_mcl_R;

R_mcl_b = zeros(size(mcl_frames,1),size(mcl_frames,2),numFrames_mcl);
R_mcl_f = zeros(size(mcl_frames,1),size(mcl_frames,2),numFrames_mcl);

for i = 1:numFrames_mcl
   R_mcl_b(:,:,i) = reshape(u_dmd_mcl_low_final(:,i), [size(mcl_frames,1),size(mcl_frames,2)]);
   R_mcl_f(:,:,i) = reshape(u_dmd_mcl_sparse_final(:,i), [size(mcl_frames,1),size(mcl_frames,2)]);
end

%%
figure(1)
imshow(mat2gray(R_mcl_b(:,:,200)));
figure(2)
imshow(mat2gray(R_mcl_f(:,:,200)));

%%
line = -189:189;
plot(zeros(numFrames_mcl,1),line,'k','Linewidth',2) % imaginary axis
hold on
plot(line,zeros(numFrames_mcl,1),'k','Linewidth',2) % real axis
plot(real(omega_mcl)*dt_mcl,imag(omega_mcl)*dt_mcl,'r.','Markersize',15)
xlabel('Re(\omega)')
ylabel('Im(\omega)')
set(gca,'FontSize',16,'Xlim',[-1.5 0.5],'Ylim',[-3 3])
%%
implay(R_mcl_b);
implay(R_mcl_f);

%% Ski Drop
M_sdl_1 = M_sdl(:,1:end-1);
M_sdl_2 = M_sdl(:,2:end);
[U_sdl,Sigma_sdl,V_sdl] = svd(M_sdl_1,'econ');
S_sdl = U_sdl'*M_sdl_2*V_sdl*diag(1./diag(Sigma_sdl));
[eV_sdl,D_sdl] = eig(S_sdl);
mu_sdl = diag(D_sdl);
omega_sdl = log(mu_sdl)/dt_sdl;
Phi_sdl = U_sdl*eV_sdl;
y0_sdl = Phi_sdl\M_sdl_1(:,1);

mode_sdl = 0;
mode_sdl_i = zeros(length(omega_sdl),1);
for i = 1:length(omega_sdl)
    if abs(real(omega_sdl(i))) < 1e-3
        mode_sdl = mode_sdl+1;
        mode_sdl_i(mode_sdl) = i;
    end
end
mode_sdl_i = mode_sdl_i(1:mode_sdl);

y0_sdl_low = y0_sdl(mode_sdl_i);
omega_sdl_low = omega_sdl(mode_sdl_i);
Phi_sdl_low = Phi_sdl(:,mode_sdl_i);

u_modes_sdl_low = zeros(length(y0_sdl_low),numFrames_sdl);
for iter = 1:numFrames_sdl
   u_modes_sdl_low(:,iter) = y0_sdl_low.*exp(omega_sdl_low*t_sdl(iter)); 
end
u_dmd_sdl_low = Phi_sdl_low*u_modes_sdl_low;

%
u_dmd_sdl_sparse = M_sdl - abs(u_dmd_sdl_low);
u_sdl_R = zeros(size(u_dmd_sdl_sparse,1),size(u_dmd_sdl_sparse,2));
for i = 1:size(u_dmd_sdl_sparse,1)
    for j = 1:size(u_dmd_sdl_sparse,2)
        if u_dmd_sdl_sparse(i,j) < 0
            u_sdl_R(i,j) = u_dmd_sdl_sparse(i,j);
        end
    end
end
u_dmd_sdl_low_final = u_sdl_R+abs(u_dmd_sdl_low);
u_dmd_sdl_sparse_final = u_dmd_sdl_sparse-u_sdl_R;

%
R_sdl_b = zeros(size(sdl_frames,1),size(sdl_frames,2),numFrames_sdl);
R_sdl_f = zeros(size(sdl_frames,1),size(sdl_frames,2),numFrames_sdl);

for i = 1:numFrames_sdl
   R_sdl_b(:,:,i) = reshape(u_dmd_sdl_low_final(:,i), [size(sdl_frames,1),size(sdl_frames,2)]);
   R_sdl_f(:,:,i) = reshape(u_dmd_sdl_sparse_final(:,i), [size(sdl_frames,1),size(sdl_frames,2)]);
end