clc; clear; close all;
addpath(genpath(pwd));

A     = load('adult.mat').X; 
b     = load('adultclass.mat').y; 
A     = Normalization(A,3);  
[d,n] = size(A);  
I     = randperm(d);
A     = A(I,:);  % shuffle samples
b     = b(I,:);  
m     = 100;     % divide samples into m groups
while 1
    idx = unique([randperm(d-2,m-1)+1 d]);
    di  = (idx-[0 idx(1:end-1)])'; 
    if min(di)>0.01*d/m; break; end
end 
 
pars.noise = 0.1;    
pars.rho   = 0.5;   
k0         = 10;
out        = FedEPM(di,n,A,b,k0,pars);
plotobj(out.objx/m) 

