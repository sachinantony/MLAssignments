clear all;
clear;
clc;
load('hm2data2.mat')
X = data;
muX = mean(X);
stdX = std(X);
repstd=repmat(stdX,97,1);
repmu=repmat(muX, 97, 1);
standardizedX = (X-repmu)./repstd;
[mean(standardizedX) ; std(standardizedX)];
Y = standardizedX(:,4);
X = standardizedX(:,1:3);
m = length(X);
eX = [ones(m,1) X];
i = 0.001;
alpha = 0.03;
iterations = 1500;
W = zeros(4,1);
WList= [];
i = 0.001;
lambda =[i];
while i< 10
    i = i*2;
    lambda = [lambda,i];
end

for i= 1 : length(lambda)
 WList = [WList,gradientDescent2(eX, Y, W, alpha, lambda(i), iterations)];
end
[Min,Index]=min(mean(WList'));
eX(:,Index)=[];