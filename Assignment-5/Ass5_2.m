clear;
clear all;
clc;

load fisheriris;

%replace species names with numerical values
for i = 1:length(species)
    if strcmp(species(i,1),'setosa')
        species(i,1) = num2cell(1);
    else 
        species(i,1) = num2cell(-1);
    end

end
C = 0.5;
numerical_species = cell2mat(species);
merged_data = horzcat(meas, numerical_species);   
clear meas;                 %delete meas
clear species;              %delete species
clear numerical_species;    %delete numerical_species
rowSize = length(merged_data(1,:));
Xi = merged_data(:,1:rowSize - 1);
Xi= zscore(Xi);
yi = merged_data(:,rowSize) ;


figure
hold on
scatter(Xi(yi==1,1),Xi(yi==1,2),'.b')
scatter(Xi(yi==-1,1),Xi(yi==-1,2),'.r')
title('Linearly seperable data')
xlabel('{x_1}'),ylabel('{x_2}')
legend('Positive class','Negative class')
hold off

n = length(Xi(:,1));
alpha=rand(1,n-1);
alpha = alpha';
prdctSum = sum(yi(1:n-1).*alpha);
alpha(n)=(-prdctSum)/yi(n);
K=Xi*Xi';
b = zeros(n,1);
for zorro = 1:100
for i = 1:n
    temp(i,:) = alpha(i)*yi(i)*Xi(i,:);
end
w=sum(temp);
for i = 1:n
    KKT1(i) = alpha(i)*(yi(i)*(w*Xi(i,:)'+b(i))-1);
end
%KKT1 = eval(KKT);
[val,i1] = max(KKT1);
alpha_old = alpha;
% for i = 1:n
%     for j = 1: n
%         tempE(i,j) = alpha(j)*yi(j)*K(i,j);
%     end
% end
for i = 1:n
    E(i) = sum(alpha.*yi.*K(i,:)')-yi(i);
end
for i = 1:n
    e(i) = E(i1)-E(i);
end
[val,i2]=max(e);
k = K(i1,i1)+K(i2,i2)-2*K(i1,i2);
alpha_new = alpha;
alpha_new(i2) = alpha_old(i2) +  yi(i2)*E(i2)/k;
if(alpha_new(i2) < 0.001)
    alpha_new(i2) = 0
end
alpha_new(i1) = alpha_old(i1) + yi(i1)*yi(i2)*(alpha_old(i2)-alpha_new(i2));
if(alpha_new(i1) < 0.001)
    alpha_new(i1) = 0
end

for i = 1:n
    if alpha(i) > 0
        b(i) = yi(i) - w*Xi(i, :)';
    end
end

bCount = 0;
for i = 1:n
    if alpha(i) ~= 0
        bCount = bCount + 1;
    end
end
totalB = sum(b) / bCount;

w = zeros(1, length(Xi(1,:)));
for i = 1:n
    w = w + alpha(i) * yi(i) * Xi(i, :);
end

% 
% b1 = b - E(i1) - yi(i1)*(alpha_new(i1)-alpha_old(i1))*K(i1,i1) - yi(i2)*(alpha_new(i2)-alpha_old(i2))*K(i2,i2);
% b2 = b - E(i2) - yi(i1)*(alpha_new(i1)-alpha_old(i1))*K(i1,i1) - yi(i2)*(alpha_new(i2)-alpha_old(i2))*K(i2,i2);
% if b1 == b2
%     b = b1;
% else
%     b = mean([b1 b2]);
% end
alpha = alpha_new
truePositiveCount = 0;
falsePositiveCount = 0;
trueNegativeCount = 0;
falseNegativeCount = 0;

for i = 1:n
    output = sign(w*Xi(i, :)' + totalB);
    if output == 1
       if output == yi(i)
           truePositiveCount = truePositiveCount + 1;
       else
           falsePositiveCount = falsePositiveCount + 1;
       end
    else
        if output == yi(i)
            trueNegativeCount = trueNegativeCount + 1;
        else
            falseNegativeCount = falseNegativeCount + 1;
        end
    end
end

truePositive = truePositiveCount / (truePositiveCount + falseNegativeCount);
falsePositive = falsePositiveCount / (falsePositiveCount + falseNegativeCount);

% break if we reached the desired accuracy (based on the param eps)
% if 1 - (truePositiveCount + trueNegativeCount) / (truePositiveCount + trueNegativeCount + falseNegativeCount + falsePositiveCount) < eps
%     break
% end
end
y=sign(w*Xi'+totalB);
Xx = w*Xi'+totalB;

plot(Xx,1);
confusionmat(yi,y)

figure 
hold on
scatter(x(y==1,1),x(y==1,2),'.b')
scatter(x(y==-1,1),x(y==-1,2),'.r')
scatter(Xsupport(Ysupport==1,1),Xsupport(Ysupport==1,2),'b')
scatter(Xsupport(Ysupport==-1,1),Xsupport(Ysupport==-1,2),'r')

syms x1 x2
fn=vpa((-bias-W*x1)/W(2),5)
fn1=vpa((-1-bias-W*x1)/W(2),5)
fn2=vpa((1-bias-W*x1)/W(2),5)
fplot(fn,'Linewidth',1);
fplot(fn1,'Linewidth',1);
fplot(fn2,'Linewidth',1);
title('Classification using SVM (Linear Kernel & SMO)')
xlabel('{x_1}'),ylabel('{x_2}')
legend('Positive class','Negative class')
hold off
