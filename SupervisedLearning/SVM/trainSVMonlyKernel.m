function [ f, alpha ] = trainSVMonlyKernel( kernel,y,C,method)
%TRAINSVM Summary of this function goes here
%   Detailed explanation goes here



weights=getWeightsMKL(kernel,y,method);
%fprintf('%f %f %f %f %f \n \n',weights(1),weights(2),weights(3),weights(4),weights(5));
K=sumUpKernels(kernel,weights);

minEig=min(eig(K));
n=size(K,1);
if (minEig<0)
    K=K-minEig*eye(n);
end

alpha=svmDual(K,y,C);

f=getHyperplaneOnlyKernel(alpha,kernel,weights,y);

end



function sumK=sumUpKernels(K,weights)

m=size(K,3);
n=size(K,1);
sumK=zeros(n,n);
for i=1:m
    sumK=sumK+weights(i)*K(:,:,i);
end

end
