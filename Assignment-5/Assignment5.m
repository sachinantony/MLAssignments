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
numerical_species = cell2mat(species);
merged_data = horzcat(meas, numerical_species);   
clear meas;                 %delete meas
clear species;              %delete species
clear numerical_species;    %delete numerical_species
rowSize = length(merged_data(1,:));
Xi = merged_data(:,1:rowSize - 1);
yi = merged_data(:,rowSize) ;
len = length(yi);
alpha=rand(1,len-1);
alpha = alpha';
prdctSum = sum(yi(1:len-1).*alpha);
alpha(len)=(-prdctSum)/yi(len);
b = 0;
w=[];
for i = 1:rowSize-1
    w=[w,alpha.*yi.*Xi(:,i)];
end
    

