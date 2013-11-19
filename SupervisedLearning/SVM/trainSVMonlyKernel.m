function [ f, alpha1 ,weights] = trainSVMonlyKernel( kernel,y,C,method)
%TRAINSVM Summary of this function goes here
%   Detailed explanation goes here


if (~strcmp(method,'multiplication'))
    weights=getWeightsMKL(kernel,y,method);
    K=sumUpKernels(kernel,weights);
else
    n=size(kernel,3);
    K=kernel(:,:,1);
    
    for i=2:n
        K=K.*kernel(:,:,i);
    end
    
    multiplication=1;

end
%weights=[1/4,1/4,1/2];
%fprintf('%f %f %f %f %f \n \n',weights(1),weights(2),weights(3),weights(4),weights(5));


minEig=min(eig(K));
n=size(K,1);
if (minEig<0)
    K=K-minEig*eye(n);
end

alpha1=svmDual(K,y,C,1);

f=getHyperplaneOnlyKernel(alpha1,kernel,weights,y);

end



function sumK=sumUpKernels(K,weights)

m=size(K,3);
n=size(K,1);
sumK=zeros(n,n);
for i=1:m
    sumK=sumK+weights(i)*K(:,:,i);
end

end
