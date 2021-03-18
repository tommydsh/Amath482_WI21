%% Clean workspace
clear all; close all; clc
format shortG

%% import files
[y1, Fs1] = audioread('GNR.m4a');
tr1 = length(y1)/Fs1; 
[y2, Fs2] = audioread('Floyd.m4a');
tr2 = length(y2)/Fs2; 
y2=y2(1:length(y2)-1);

%% Part 1
% GNR
n1 = length(y1);
t1 = linspace(0,tr1,n1+1); t1 = t1(1:n1);
k1 = (1/tr1)*[0:n1/2-1 -n1/2:-1];
ks1 = fftshift(k1);
a = 100;
tau1 = 0:0.1:tr1;
gnr_notes=zeros(length(tau1),1);
for j = 1:length(tau1)
   g1 = exp(-a*(t1 - tau1(j)).^2);
   yg1 = g1.'.*y1;
   yg1t = fft(yg1);
   [M,i]=max(abs(yg1t));
   gnr_notes(j)=k1(i);
   yg1t_spec(:,j) = fftshift(abs(yg1t));
end
round(gnr_notes,2)
figure(1)
pcolor(tau1,ks1,yg1t_spec)
shading interp
set(gca,'ylim',[0 1500],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency')

%%
% Floyd
n2 = length(y2);
t2 = linspace(0,tr2,n2+1); t2 = t2(1:n2);
k2 = (1/tr2)*[0:n2/2-1 -n2/2:-1];
ks2 = fftshift(k2); 
a = 100;
tau2 = 0:1:tr2;
gau_tau2=0.1;
floyd_notes=zeros(length(tau2),1);
for j = 1:length(tau2)
   g2 = exp(-a*(t2 - tau2(j)).^2); % Window function
   yg2 = g2.'.*y2;
   yg2t = fft(yg2);
   yg2t_spec(:,j) = fftshift(abs(yg2t)); % We don't want to scale it
   
   [M,idx2]=max(abs(yg2t));
   floyd_notes(j)=k2(idx2);
   
   gau_filter=exp(-gau_tau2*(k2-abs(k2(idx2))).^2);
   yg2t_filtered=gau_filter'.*yg2t;
   yg2t_spec_filtered(:,j)=fftshift(abs(yg2t_filtered));
end
% k2(i)
figure(2)
pcolor(tau2,ks2,yg2t_spec)
shading interp
set(gca,'ylim',[0 300],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency')
figure(3)
pcolor(tau2,ks2,yg2t_spec_filtered)
shading interp
set(gca,'ylim',[0 300],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency')
floyd_notes

%% Part 3
% Floyd
n3 = length(y2);
t3 = linspace(0,tr2,n3+1); t3 = t3(1:n3);
k3 = (1/tr2)*[0:n3/2-1 -n3/2:-1];
ks3 = fftshift(k3); 
a = 100;
tau3 = 0:1:tr2;
freq_bass=[82.407;92.499;97.999;110.0;123.47];
floyd_notes_guitar=zeros(length(tau3),1);
gau_tau3=0.1
for j = 1:length(tau3)
   g3 = exp(-a*(t3 - tau3(j)).^2); 
   yg3 = g3.'.*y2;
   yg3t = fft(yg3);
   
   for i=1:length(yg3t)
       for k=1:5
           if abs(abs(k3(i))-freq_bass(k))<=5
               yg3t(i)=0;
               break;
           end
       end
   end
%    yg3t_spec(:,j) = fftshift(abs(yg3t));
   [M,idx3]=max(abs(yg3t));
   floyd_notes_guitar(j)=abs(k3(idx3));
   
   gau_filter=exp(-gau_tau3*(k3-abs(k3(idx3))).^2);
   yg3t_filtered=gau_filter'.*yg3t;
   yg3t_spec_filtered(:,j)=fftshift(abs(yg3t_filtered));
end
floyd_notes_guitar
figure(3)
pcolor(tau3,ks3,yg3t_spec_filtered)
shading interp
set(gca,'ylim',[0 1000],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency')