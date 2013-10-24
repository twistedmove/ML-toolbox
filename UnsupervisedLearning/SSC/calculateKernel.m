function [ kernelMatrix] = calculateKernel( X, kernel )
%CALCULATEKERNEL calculates 3 dim kernel matrix from multiple data and
%   kernels.
%   
%   Input:
%       X               -       m x n cell array with the data; n denotes
%           number of training samples
%       kernel          -       m x 1 cell array of function handles
%           representing kernels
%
%
%   Output:
%
%       K               -       (unnormalized)positive definite kernel 
%                                   matrix which is a convex combination of
%                                   the given kernels
%
%   author: Ivan Bogun
%   data  : June 29, 2013

% Note : won't work for ordinary arrays



[~,n]=size(X);
m=length(kernel);
kernelMatrix=zeros(n,n,m);

for k=1:m
    
    if (isa(kernel{k},'double'))
       kernelMatrix(:,:,k)=kernel{k};
       kernelMatrix(:,:,k)=normalizeKernel(kernelMatrix(:,:,k));
       continue;
       
    end
    
    for i=1:n
        for j=1:n
            if (j>=i)
                if (iscell(kernel)==1)
            kernelMatrix(i,j,k)=kernel{k}(X{k,i},X{k,j});
                else
                    kernelMatrix(i,j,k)=kernel(X{k,i},X{k,j});
                end
            else
              kernelMatrix(i,j,k)=kernelMatrix(j,i,k);
            end
        end
    end
    
    kernelMatrix(:,:,k)=normalizeKernel(kernelMatrix(:,:,k));
end



end


function K=normalizeKernel(K)

n=size(K,1);

for i=1:n
    for j=1:n
        K(i,j)=K(i,j)/sqrt(K(i,i)*K(j,j));
    end
end

end
