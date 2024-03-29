function [ hyperplane ] = getHyperplane( X,alpha,kernelFn,weights,y)
%GETHYPERPLANE From solution of the dual svm get the function representing
%separating hyperplane.
%
%   Input:
%       X           -           1 x m cell array with m data points
%       alpha       -           solution of the dual soft-margin SVM
%       kernelFn    -           1 x P cell array with kernel functions
%       weights     -           1 x P array of kernel weights
%       y           -           m x 1 label matrix
%
%   Output:
%       hyperplane  -           function handle which represents separating
%           hyperplane
%
%   author: Ivan Bogun
%   date  : June 29, 2013
%
% Note : works only with cell arrays
%

f=@(x) sumOverSupportVectors(X,alpha,kernelFn,weights,y,x);

if nargin<5
    weights=1;
end

n=length(y);

maxSupport=0;
minSupport=0;

for j=1:n
    if (y(j)==-1)
        maxSupport=max(maxSupport,f(X(:,j)));
    else
        minSupport=min(minSupport,f(X(:,j)));
    end
end

b=-(maxSupport+minSupport)/2;

hyperplane=@(x) f(x)+b;

    function res=sumOverSupportVectors(X,alpha,kernelFn,weights,y,x)
        idx=find(alpha~=0);
        
        n=length(idx);
        m=length(kernelFn);
        res=0;
        
        for i=1:n
            ker=0;
            for j=1:m
                if (isa(kernelFn{j},'double'))
                    ker=ker+weights(j)*x{j}(idx(i));
                else
                    ker=ker+weights(j)*kernelFn{j}(X{j,idx(i)},x{j});
                end
                
                
            end
            res=res+alpha(idx(i))*y(idx(i))*ker;
        end
        
        
    end


end

