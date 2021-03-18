%% HW4
% Clean workspace
clear all; close all; clc

%% Import Data
[images, labels] = mnist_parse('train-images-idx3-ubyte', 'train-labels-idx1-ubyte');
[test_images, test_labels] = mnist_parse('t10k-images-idx3-ubyte', 't10k-labels-idx1-ubyte');

M = zeros(784,60000);
for i=1:60000
   M(:,i)=reshape(images(:,:,i),[784,1]); 
end
Mtest = zeros(784,10000);
for i=1:10000
   Mtest(:,i)=reshape(test_images(:,:,i),[784,1]); 
end

[m,n]=size(M);
mn=mean(M,2);
M=M-repmat(mn,1,n);
[U,S,V] = svd(M, 'econ');

sigma = diag(S);
energy = zeros(1,length(sigma));
for j = 1:length(sigma)
    energy(j) = sigma(j)^2/sum(sigma.^2);
end
figure(1)
for j=1:10
   subplot(2,5,j)
   ut1 = reshape(U(:,j),28,28);
   ut2 = rescale(ut1);
   imshow(ut2)
end

% Singular Value Spectrum
figure(2)
plot(energy,'ko','Linewidth',2)
title('Energy')

% Projection onto 3 V-modes
figure(3)
for label=0:9
    label_indices = find(labels == label);
    plot3(V(label_indices, 4)*100, V(label_indices, 5)*100, V(label_indices, 6)*100,...
        'o', 'DisplayName', sprintf('%i',label), 'Linewidth', 2)
    hold on
end
xlabel('2nd V-Mode'), ylabel('3rd V-Mode'), zlabel('5th V-Mode')
title('Projection onto V-modes 2, 3, 5')
legend
set(gca,'Fontsize', 14)
clc

%% 2 digits
LDAsuc = zeros(10);
LDAsuct = zeros(10);
for i = 1:10
   train_data{i} = M(:,find(labels==i-1));
   test_data{i} = Mtest(:,find(test_labels==i-1));
end

%%
feature=10;
for i = 1:9
    for j = i+1:10
        num1 = i-1;
        num2 = j-1;
        num1_train = train_data{num1+1};
        num2_train = train_data{num2+1};
        num1_test = test_data{num1+1};
        num2_test = test_data{num2+1};

        % LDA
        [Ua,Sa,Va,threshold,w,sort1,sort2] = dc_trainer(num1_train,num2_train,feature);

        num1_len=size(num1_train,2);
        num2_len=size(num2_train,2);
        num1_tlen=size(num1_test,2);
        num2_tlen=size(num2_test,2);

        TrainM = [num1_train num2_train];
        TrainA = [zeros(1,num1_len) ones(1,num2_len)];
        TestM = [num1_test num2_test];
        TestA = [zeros(1,num1_tlen) ones(1,num2_tlen)];
        TrainMat = Ua'*TrainM;
        TestMat = Ua'*TestM;
        pval = w'*TrainMat;
        pvalt = w'*TestMat;

        ResVec = (pval > threshold);
        ResVec_test = (pvalt > threshold);
        err = abs(ResVec - TrainA);
        errt = abs(ResVec_test - TestA);
        errNum = sum(err);
        errtNum = sum(errt);
        LDAsuc(i,j) = 1 - errNum/(num1_len+num2_len);
        LDAsuct(i,j) = 1 - errtNum/(num1_tlen+num2_tlen);
    end
end

%%
% Decision Tree
num1 = 7;
num2 = 9;
num1_train = train_data{num1+1};
num2_train = train_data{num2+1};
num1_test = test_data{num1+1};
num2_test = test_data{num2+1};

train_data = [num1_train num2_train];
result1 = ones(1,size(num1_train,2)).*num1;
result2 = ones(1,size(num2_train,2)).*num2;
train_result = [result1 result2];
test_data = [num1_test num2_test];
tresult1 = ones(1,size(num1_test,2)).*num1;
tresult2 = ones(1,size(num2_test,2)).*num2;
test_result = [tresult1 tresult2];

tree=fitctree(train_data',train_result,'CrossVal','on');
% view(tree.Trained{1},'Mode','graph');
tree_train = predict(tree.Trained{1},train_data');
tree_test = predict(tree.Trained{1},test_data');
dterr = 0;
dterrt = 0;
for i = 1:length(train_result)
   if train_result(i) ~= tree_train(i)
      dterr = dterr+1; 
   end
end
for i = 1:length(test_result)
   if test_result(i) ~= tree_test(i)
      dterrt = dterrt+1; 
   end
end
DTsuc = 1-dterr/length(train_result)
DTsuct = 1-dterrt/length(test_result)
classError = kfoldLoss(tree);

% SVM
Mdl = fitcsvm(train_data',train_result);
train_labels_predict = predict(Mdl,train_data');
test_labels_predict = predict(Mdl,test_data');
svmerr = 0;
svmerrt = 0;
for i=1:length(train_result)
    if train_result(i) ~= train_labels_predict(i)
       svmerr = svmerr+1; 
    end
end
for i=1:length(test_result)
    if test_result(i) ~= test_labels_predict(i)
       svmerrt = svmerrt+1; 
    end
end
SVMsuc = 1-svmerr/length(train_result)
SVMsuct = 1-svmerrt/length(test_result)

%% 3 digits LDA

num1 = 0;
num2 = 1;
num3 = 4;
num1_train = train_data{num1+1};
num2_train = train_data{num2+1};
num3_train = train_data{num3+1};
num1_test = test_data{num1+1};
num2_test = test_data{num2+1};
num3_test = test_data{num3+1};

num1_wave = dc_wavelet(num1_train);
num2_wave = dc_wavelet(num2_train);
num3_wave = dc_wavelet(num3_train);
feature = 10;

[U3,S3,V3,threshold_low,threshold_high,w3,mini,midi,maxi] = dc_trainer3(num1_train,num2_train,num3_train,feature);
num1_len=size(num1_test,2);
num2_len=size(num2_test,2);
num3_len=size(num3_test,2);
position = [mini,midi,maxi]
for i=1:3
    newp(i) = find(position==i);
end
newp = newp-1;

TestM3 = [num1_test num2_test num3_test];
TestA3 = [ones(1,num1_len).*newp(1) ones(1,num2_len).*newp(2) ones(1,num3_len).*newp(3)];
TestMat3 = U3'*TestM3; 
pval = w3'*TestMat3;

ResVec = zeros(1,num1_len+num2_len+num3_len);
yes=0;
for i=1:num1_len+num2_len+num3_len
    if pval(i) < threshold_low
        ResVec(i) = 0;
    elseif pval(i) > threshold_high
        ResVec(i) = 2;
    else
        ResVec(i) = 1;
    end
    if ResVec(i) == TestA3(i)
        yes=yes+1;
    end
end
sucRate = yes/(num1_len+num2_len+num3_len)
