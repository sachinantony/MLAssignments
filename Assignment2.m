clear all;
clear;
clc;

%loading the input matrix
load('hm2data2.mat')
X = data;
Y = X(:,4);
%calculating the mean and Standard Deviation
muX = mean(X);
stdX = std(X);
m=length(X)
repstd=repmat(stdX,m,1);
repmu=repmat(muX, m, 1);
standardizedX = (X-repmu)./repstd;
[mean(standardizedX) ; std(standardizedX)];

%Separating  and Y
%Y = standardizedX(:,4);
X = standardizedX(:,1:3);
m = length(X);
eX = [ones(m,1) X];
alpha = 0.01;
iterations = 1500;
W = zeros(4,1);
halfSize = ceil(m/2);
trDataX = eX(1:halfSize,:);
trDataY = Y(1:halfSize,:);
tsDataX = eX(halfSize+1:m,:);
tsDataY = Y(halfSize+1:m,:);




%------------------------Regularization----------------------------------
WList= [];
i = 0.01;
% lambda = [i];
lambda =logspace(-2,1,25);
% while i< 10
%     i = i*2;
%     lambda = [lambda,i];
% end
%lambda = logspace()
for i= 1 : length(lambda)
 WList = [WList,gradientDescent2(trDataX, trDataY, W, alpha, lambda(i), iterations)];
end
WList = WList';

%-----------------Get lambda-------------------
costList = [];
minCost = computeCost2(trDataX, trDataY,WList(1,:)',lambda(1))
for i = 1:length(lambda)
    costList = [costList, computeCost2(trDataX, trDataY,WList(i,:)',lambda(i))];
end
costList
[minValue,minIndex] = min(costList);
minIndex
lambdaVal = lambda(minIndex)
% minIndex=1;
% for i = 1:length(WList)
%     omega = WList(i,:);
%     omega = omega';
%     J = computeCost2(trDataX, trDataY,omega,lambda(i));
%     if computeCost2(trDataX, trDataY, omega,lambda(i)) < minCost
%         minCost = J;
%         minIndex = i;
%     end
% end
% currentLambda = lambda(minIndex);


[Min,Index]=min(mean(WList));
trDataX;
trDataX(:,Index)=[];
eX(:,Index)=[];
trainLength = ceil(m/2);
testLength = m - trainLength;

%-----------------Cross Validation Starts------------------------------
W = zeros(3,1);
% subsetSize = round(m/3);
% XSet1 = eX(1:subsetSize,:);
% XSet2 = eX(subsetSize+1:subsetSize*2,:);
% XSet3 = eX(2*subsetSize+1:m,:);
% YSet1 = Y(1:subsetSize,:);
% YSet2 = Y(subsetSize+1:subsetSize*2,:);
% YSet3 = Y(2*subsetSize+1:m,:);
% 
% Set1W = gradientDescentB(XSet1, YSet1, W, alpha, iterations);
% Set2W = gradientDescentB(XSet2, YSet2, W, alpha, iterations);
% Set3W = gradientDescentB(XSet3, YSet3, W, alpha, iterations);
% 
% costSet1 = [computeCostB(XSet2, YSet2, Set1W);computeCostB(XSet3, YSet3, Set1W)]
% costSet2 = [computeCostB(XSet1, YSet1, Set2W);computeCostB(XSet3, YSet3, Set1W)]
% costSet3 = [computeCostB(XSet1, YSet1, Set3W);computeCostB(XSet2, YSet2, Set1W)]
% 
% mean([costSet1;costSet2;costSet3;])

%cvsvm = fitcecoc(trDataX,trDataY,'Kfold',3);
%Yhat = kfoldPredict(cvsvm);
deg1Cost = crossValidation3Fold(eX, Y, W, alpha, iterations);
deg2Polynomial = [eX,eX(:,2:end).^2];
W=zeros(5,1);
deg2Cost = crossValidation3Fold(deg2Polynomial, Y, W, alpha, iterations);
if(deg1Cost > deg2Cost)
    eX = deg2Polynomial;
end
modelingError=[];
genError = [];
iList = []
for i = 1:100
    iList = [iList,i];
    currentX = [];
    currentY = [];
    randList = randperm(m);
    for j = 1:m
        currentX = [currentX;eX(randList(j),:)];
        currentY = [currentY;Y(randList(j),:)];
    end

    trSetX = currentX(1:halfSize,:);
    trSetY = currentY(1:halfSize,:);
    tstSetX = currentX(halfSize:end,:);
    tstSetY = currentY(halfSize:end,:);
    currentW = gradientDescentB(trSetX,trSetY,W,alpha,iterations);
    modelingError = [modelingError,computeCost2(trSetX,trSetY,currentW,lambdaVal)];
    genError = [genError,computeCost2(tstSetX,tstSetY,currentW,lambdaVal)];
end
figure
plot(iList,modelingError);
hold on;
plot(iList,genError)
hold off;
    
    






