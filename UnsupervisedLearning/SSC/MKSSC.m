function [predicted] = MKSSC( X,n,K,lambda,kernels,kernelWeights)
%SSC Sparse subspace clustering algorithm
%   Multiple kernel sparse subspace clustering algorithm
%
% Input:
%   X               -       mxn data matrix
%   groundTruth     -       nx1 ground truth
%   n               -       how many subspaces to seek
%   K               -       number of coefficients to take for the
%       projection
%   lambda          -       regularization parameter
%   kernel          -       function handle to the kernel, if using kernel
%       SSC algorithm
%
%
%  Output:
%   predicted       -       predicted labels
%   confusionMatrix -       confusion matrix of the clustering
%   Missrate        -       missclassification rate of the clustering
%
%   author: Ivan Bogun
%   date  : June 10, 2013
%
%  credit: code is adapted from the one by Ehsan Elhamifar



m=length(kernels);
if (nargin<6)
    kernelWeights=ones(1,m)*(1/m);
end

allKernel=calculateKernel(X,kernels);
combinedKernel=kernelSum(allKernel,kernelWeights);

for i=1:10
    combinedKernel=kernelSum(allKernel,kernelWeights);
    CMat = calculateSparseCoefficients(X,lambda,combinedKernel);
    CKSym = BuildAdjacency(CMat,K);
    [predicted , ~, ~] = SpectralClustering(CKSym,n);
            
    load('groundTruth');
    [Missrate, confusionMatrix,predicted] = Misclassification(predicted,groundTruth);
    display(min(Missrate));
    
    
    kernelWeights=solveForKernelWeights(CMat,allKernel);
    display(kernelWeights);
end

end


function [ K]=kernelSum(kernelMatrix,kernelWeights)

if (any(kernelWeights<0)|| sum(kernelWeights)~=1)
    fprintf('kernel weights have to be positive and should sum up to one');
    return;
end
m=size(kernelMatrix,3);
n=size(kernelMatrix,1);

K=zeros(n,n);

for k=1:m
    K=K+kernelWeights(k)*kernelMatrix(:,:,k);
end

minEig=min(eig(K));

if (minEig<0)
    K=K-eye(n)*minEig;
end

end


function [gamma]=solveForKernelWeights(CMat,kernelMatrix)

n=size(CMat);
m=size(kernelMatrix,3);

% calculate K1,K2,K3
K1=zeros(n,m);
K2=zeros(n,m);
K3=zeros(n,m);


for j=1:m
    for i=1:n
        
        c_i=CMat(:,i);
        
        K1(i,j)=kernelMatrix(i,i,j);
        
        K2(i,j)=-2*c_i'*kernelMatrix(:,i,j);
        
        K3(i,j)=quad_form(c_i,kernelMatrix(:,:,j));
    end
end


K=K1+K2+K3;

cvx_begin
cvx_quiet true
cvx_precision high
variable gamma(m)
minimize ( sum(K*gamma))
subject to
    sum(gamma)==1;
    gamma>=0;
cvx_end




end
