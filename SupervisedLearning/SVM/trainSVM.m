function [ f, alpha ] = trainSVM( X,y,kernel,C,method)
%TRAINSVM Summary of this function goes here
%   Detailed explanation goes here


K=calculateKernel(X,kernel);

if (nargin>4)
   weights=getWeightsMKL(K,y,method);
   K=sumUpKernels(K,weights);
   
   %kernel=@(x) kernelValue(kernel,weights,x);
   
end

minEig=min(eig(K));
n=size(K,1);
if (minEig<0)
K=K-minEig*eye(n);
end

alpha=svmDual(K,y,C);

f=getHyperplane(X,alpha,kernel,weights,y);

end



function sumK=sumUpKernels(K,weights)

m=size(K,3);
n=size(K,1);
sumK=zeros(n,n);
for i=1:m
    sumK=sumK+weights(i)*K(:,:,i);
end

end

function res=kernelValue(kernel,weights,x)

res=0;
m=length(weights);

for i=1:m
   res=res+weights(i)*kernel{i}(x); 
end

end