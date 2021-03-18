%% Project 1
% Clean workspace
clear all; close all; clc
format short
load subdata.mat % Imports the data as the 262144x49 (space by time) matrix called subdata
L = 10; % spatial domain
n = 64; % Fourier modes
x2 = linspace(-L,L,n+1); x = x2(1:n); y = x; z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

%% Avergaing Frequencies
Unave = zeros(64,64,64);
for j=1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n);
    Unfftn=fftn(Un);
    Unave = Unave + Unfftn;
end
Unave = abs(fftshift(Unave))/49;
[M,idx]=max(abs(Unave),[],'all','linear');
[y0,x0,z0]=ind2sub(size(Unave),idx);
v=[x0,y0,z0];
v=v.*(20/64)-10
figure(1)
isosurface(Kx,Ky,Kz,abs(Unave)/M,0.7)
axis([-20 20 -20 20 -20 20]), grid on, drawnow,
xlabel('kx'),ylabel('ky'),zlabel('kz')

%% Filter around center frequency to find the path of submarine
tau=0.2;
coor=zeros(49,3);
filter = exp(-tau*((Kx-v(1)).^2 +(Ky-v(2)).^2+(Kz-v(3)).^2));
% filter = exp(-tau*((Kx-5).^2 +(Ky+7).^2+(Kz-2).^2));
for k=1:49
    Un(:,:,:)=reshape(subdata(:,k),n,n,n);
    Unfftn=fftn(Un);
    Unfftn=filter.*fftshift(Unfftn);
    Unfftn=ifftshift(Unfftn);
    Unf=ifftn(Unfftn);
    [M,idxn] = max(abs(Unf),[],'all','linear');
    [b,a,c]=ind2sub(size(Unf),idxn);
    coor(k,1)=a*20/64-10;
    coor(k,2)=b*20/64-10;
    coor(k,3)=c*20/64-10;
end
figure(2)
plot3(coor(:,1),coor(:,2),coor(:,3),'bo-','Linewidth',1)
axis([-10 10 -10 10 -10 10]), grid on, drawnow
xlabel('x'),ylabel('y'),zlabel('z')
%% x,y coordinates for the submarine
coor(:,[1,2])
%%
% randn('seed',314);
x=[-2,-1,0,1,2];
y=[2,-10,0,2,1];
p = polyfit(x,y,2)
 
y1 = polyval(p,x);
figure
plot(x,y,'o')
hold on
plot(x,y1)
hold off
