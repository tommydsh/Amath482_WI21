%% Clean workspace
clear all; close all; clc
format shortG

%%
load('cam1_1.mat')
% implay(vidFrames1_1)
load('cam2_1.mat')
% implay(vidFrames2_1))
load('cam3_1.mat')
% implay(vidFrames3_1)

%%
% 1_1
numFrames1_1 = size(vidFrames1_1,4);
x1_1 = zeros(1,numFrames1_1);
y1_1 = zeros(1,numFrames1_1);
for j = 1:numFrames1_1
    X1_1 = rgb2gray(vidFrames1_1(:,:,:,j));
    X1_1=im2double(X1_1);
    X1_1(:,1:300)=0;
    X1_1(:,400:end)=0;
    X1_1(1:200,:)=0;
    [M,i]=max(X1_1(:));
    [y0,x0]=ind2sub(size(X1_1),i);
    x1_1(j)=x0;
    y1_1(j)=y0;
%    imshow(X1_1);drawnow
end

%2_1
numFrames2_1 = size(vidFrames2_1,4);
x2_1 = zeros(1,numFrames2_1);
y2_1 = zeros(1,numFrames2_1);
for j = 1:numFrames2_1
    X2_1 = rgb2gray(vidFrames2_1(:,:,:,j));
    X2_1=im2double(X2_1);
    X2_1(:,1:250)=0;
    X2_1(:,350:end)=0;
    X2_1(1:80,:)=0;
    X2_1(400:end,:)=0;
    [M,i]=max(X2_1(:));
    [y0,x0]=ind2sub(size(X2_1),i);
    x2_1(j)=x0;
    y2_1(j)=y0;
%    imshow(X2_1);drawnow
end

% 3_1
numFrames3_1 = size(vidFrames3_1,4);
x3_1 = zeros(1,numFrames3_1);
y3_1 = zeros(1,numFrames3_1);
for j = 1:numFrames3_1
    X3_1 = rgb2gray(vidFrames3_1(:,:,:,j));
    X3_1=im2double(X3_1);
    X3_1(1:240,:)=0;
    X3_1(330:end,:)=0;
    X3_1(:,1:250)=0;
    X3_1(:,480:end)=0;
    [M,i]=max(X3_1(:));
    [y0,x0]=ind2sub(size(X3_1),i);
    x3_1(j)=x0;
    y3_1(j)=y0;
%     imshow(X3_1);drawnow
end

x2_1=x2_1(1:226);
y2_1=y2_1(1:226);
x3_1=x3_1(1:226);
y3_1=y3_1(1:226);

A1=[x1_1;y1_1;x2_1;y2_1;x3_1;y3_1];
[m,n]=size(A1); % compute data size
mn=mean(A1,2); % compute mean for each row
A1=A1-repmat(mn,1,n); 
A11=A1/sqrt(n-1);
[U1,S1,V1]=svd(A11,'econ');
Y1=U1'*A1;

for i=1:6
    subplot(3,2,i)
    plot(A1(i,:))
end

figure(2)
plot(A1(1,:)), hold on
plot(A1(4,:))
title('Ideal Case')
xlabel('Time Frame')
ylabel('displacement')
legend('horizontal','vertical')

%%
load('cam1_2.mat')
% implay(vidFrames1_2)
load('cam2_2.mat')
% implay(vidFrames2_2)
load('cam3_2.mat')
% implay(vidFrames3_2)

%%
% 1_2
numFrames1_2 = size(vidFrames1_2,4);
x1_2 = zeros(1,numFrames1_2);
y1_2 = zeros(1,numFrames1_2);
for j = 1:numFrames1_2
    X1_2 = rgb2gray(vidFrames1_2(:,:,:,j));
    X1_2(:,1:300)=0;
    X1_2(:,400:end)=0;
    X1_2(1:200,:)=0;
    [M,i]=max(abs(X1_2),[],'all','linear');
    [y0,x0]=ind2sub(size(X1_2),i);
    x1_2(j)=x0;
    y1_2(j)=y0;
%    imshow(X1_2);drawnow
end

% 2_2
numFrames2_2 = size(vidFrames2_2,4);
x2_2 = zeros(1,numFrames2_2);
y2_2 = zeros(1,numFrames2_2);
for j = 1:numFrames2_2
    X2_2 = rgb2gray(vidFrames2_2(:,:,:,j));
    X2_2(:,1:200)=0;
    X2_2(:,400:end)=0;
    X2_2(1:50,:)=0;
    X2_2(350:end,:)=0;
    [M,i]=max(abs(X2_2),[],'all','linear');
    [y0,x0]=ind2sub(size(X2_2),i);
    x2_2(j)=x0;
    y2_2(j)=y0;
%    imshow(X2_2);drawnow
end

% 3_2
numFrames3_2 = size(vidFrames3_2,4);
x3_2 = zeros(1,numFrames3_2);
y3_2 = zeros(1,numFrames3_2);
for j = 1:numFrames3_2
    X3_2 = rgb2gray(vidFrames3_2(:,:,:,j));
    X3_2(1:150,:)=0;
    X3_2(300:end,:)=0;
    X3_2(:,1:230)=0;
    X3_2(:,500:end)=0;
    [M,i]=max(abs(X3_2),[],'all','linear');
    [y0,x0]=ind2sub(size(X3_2),i);
    x3_2(j)=x0;
    y3_2(j)=y0;
%     imshow(X3_2);drawnow
end

x2_2=x2_2(1:314);
y2_2=y2_2(1:314);
x3_2=x3_2(1:314);
y3_2=y3_2(1:314);

A2=[x1_2;y1_2;x2_2;y2_2;x3_2;y3_2];

[m,n]=size(A2); % compute data size
mn=mean(A2,2); % compute mean for each row
A2=A2-repmat(mn,1,n); % subtract mean
Cx2=(1/(n-1))*A2*A2'; 
A22=Cx2/sqrt(n-1);
[U2,S2,V2]=svd(A22,'econ');
Y2=U2'*A2;

for i=1:6
    subplot(3,2,i)
    plot(A2(i,:))
end

figure(2)
plot(A2(1,:)),hold on
plot(A2(2,:))
title('Noisy Case')
xlabel('Time Frame')
ylabel('displacement')
legend('horizontal','vertical')
legend('PC1','PC2')
%%
load('cam1_3.mat')
% implay(vidFrames1_3)
load('cam2_3.mat')
% implay(vidFrames2_3)
load('cam3_3.mat')
% implay(vidFrames3_3)

%%
% 1_3
numFrames1_3 = size(vidFrames1_3,4);
x1_3 = zeros(1,numFrames1_3);
y1_3 = zeros(1,numFrames1_3);
for j = 1:numFrames1_3
    X1_3 = rgb2gray(vidFrames1_3(:,:,:,j));
    X1_3(:,1:300)=0;
    X1_3(:,400:end)=0;
    X1_3(1:200,:)=0;
    [M,i]=max(abs(X1_3),[],'all','linear');
    [y0,x0]=ind2sub(size(X1_3),i);
    x1_3(j)=x0;
    y1_3(j)=y0;
%    imshow(X1_3);drawnow
end

% 2_3
numFrames2_3 = size(vidFrames2_3,4);
x2_3 = zeros(1,numFrames2_3);
y2_3 = zeros(1,numFrames2_3);
for j = 1:numFrames2_3
    X2_3 = rgb2gray(vidFrames2_3(:,:,:,j));
    X2_3(:,1:200)=0;
    X2_3(:,400:end)=0;
    X2_3(1:50,:)=0;
    X2_3(350:end,:)=0;
    [M,i]=max(abs(X2_3),[],'all','linear');
    [y0,x0]=ind2sub(size(X2_3),i);
    x2_3(j)=x0;
    y2_3(j)=y0;
%    imshow(X2_3);drawnow
end

% 3_3
numFrames3_3 = size(vidFrames3_3,4);
x3_3 = zeros(1,numFrames3_3);
y3_3 = zeros(1,numFrames3_3);
for j = 1:numFrames3_3
    X3_3 = rgb2gray(vidFrames3_3(:,:,:,j));
    X3_3(1:150,:)=0;
    X3_3(300:end,:)=0;
    X3_3(:,1:230)=0;
    X3_3(:,500:end)=0;
    [M,i]=max(abs(X3_3),[],'all','linear');
    [y0,x0]=ind2sub(size(X3_3),i);
    x3_3(j)=x0;
    y3_3(j)=y0;
%     imshow(X3_3);drawnow
end

x1_3=x1_3(1:237);
y1_3=y1_3(1:237);
x2_3=x2_3(1:237);
y2_3=y2_3(1:237);
x3_3=x3_3(1:237);
y3_3=y3_3(1:237);

A3=[x1_3;y1_3;x2_3;y2_3;x3_3;y3_3];

[m,n]=size(A3); % compute data size
mn=mean(A3,2); % compute mean for each row
A3=A3-repmat(mn,1,n); % subtract mean
Cx3=(1/(n-1))*A3*A3'; 
A33=Cx3/sqrt(n-1);
[U3,S3,V3]=svd(A33,'econ');
Y3=U3'*A3;

for i=1:6
    subplot(3,2,i)
    plot(A3(i,:))
end

figure(2)
plot(A3(1,:)),hold on
plot(A3(2,:))
title('Horizontal Displacement')
xlabel('Time Frame')
ylabel('displacement')
legend('horizontal','vertical')
legend('PC1','PC2')

%% 
load('cam1_4.mat')
% implay(vidFrames1_4)
load('cam2_4.mat')
% implay(vidFrames2_4)
load('cam3_4.mat')
% implay(vidFrames3_4)

%%
% 1_4
numFrames1_4 = size(vidFrames1_4,4);
x1_4 = zeros(1,numFrames1_4);
y1_4 = zeros(1,numFrames1_4);
for j = 1:numFrames1_4
    X1_4 = rgb2gray(vidFrames1_4(:,:,:,j));
    X1_4(:,1:300)=0;
    X1_4(:,400:end)=0;
    X1_4(1:200,:)=0;
    [M,i]=max(abs(X1_4),[],'all','linear');
    [y0,x0]=ind2sub(size(X1_4),i);
    x1_4(j)=x0;
    y1_4(j)=y0;
%    imshow(X1_3);drawnow
end

% 2_4
numFrames2_4 = size(vidFrames2_4,4);
x2_4 = zeros(1,numFrames2_4);
y2_4 = zeros(1,numFrames2_4);
for j = 1:numFrames2_4
    X2_4 = rgb2gray(vidFrames2_4(:,:,:,j));
    X2_4(:,1:200)=0;
    X2_4(:,400:end)=0;
    X2_4(1:50,:)=0;
    X2_4(350:end,:)=0;
    [M,i]=max(abs(X2_4),[],'all','linear');
    [y0,x0]=ind2sub(size(X2_4),i);
    x2_4(j)=x0;
    y2_4(j)=y0;
%    imshow(X2_4);drawnow
end

% 3_4
numFrames3_4 = size(vidFrames3_4,4);
x3_4 = zeros(1,numFrames3_4);
y3_4 = zeros(1,numFrames3_4);
for j = 1:numFrames3_4
    X3_4 = rgb2gray(vidFrames3_4(:,:,:,j));
    X3_4(1:150,:)=0;
    X3_4(300:end,:)=0;
    X3_4(:,1:230)=0;
    X3_4(:,500:end)=0;
    [M,i]=max(abs(X3_4),[],'all','linear');
    [y0,x0]=ind2sub(size(X3_4),i);
    x3_4(j)=x0;
    y3_4(j)=y0;
%     imshow(X3_4);drawnow
end

x1_4=x1_4(1:237);
y1_4=y1_4(1:237);
x2_4=x2_4(1:237);
y2_4=y2_4(1:237);
x3_4=x3_4(1:237);
y3_4=y3_4(1:237);

A4=[x1_4;y1_4;x2_4;y2_4;x3_4;y3_4];

[m,n]=size(A4); % compute data size
mn=mean(A4,2); % compute mean for each row
A4=A4-repmat(mn,1,n); % subtract mean
Cx4=(1/(n-1))*A4*A4'; 
A44=Cx4/sqrt(n-1);
[U4,S4,V4]=svd(A44,'econ');
Y4=U4'*A4;

for i=1:6
    subplot(3,2,i)
    plot(A4(i,:))
end

figure(2)
plot(A4(1,:)),hold on
plot(A4(2,:))
title('Horizontal Displacement and Rotation')
xlabel('Time Frame')
ylabel('displacement')
legend('horizontal','vertical')
legend('PC1','PC2')