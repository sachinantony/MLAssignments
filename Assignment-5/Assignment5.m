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


n = 50; 
rand ( 'seed' ,2) ; 
Xi = 6* rand (n ,2) ;
q = 0;
Xi = [Xi, 4* rand(n,q) ];
[n,p] = size (Xi) ;
bt = -6; 
wt = [4 ; -1]; 
yi = sign (wt (1) * Xi (: ,1) + wt (2) * Xi (: ,2) + bt); 
len = length(yi);
alpha=randperm(len-1);
alpha = alpha';
prdctSum = sum(yi(1:len-1).*alpha);
alpha(len)=(-prdctSum)/yi(len);

