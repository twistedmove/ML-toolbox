function [predicted, confusionMatrix, Missrate] = SSC( X,groundTruth,n,K,lambda,kernel)
%SSC Sparse subspace clustering algorithm
%   Detailed explanation goes here

% Input:
%   X               -       DxN data matrix
%   groundTruth     -       Nx1 ground truth
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

if nargin==6
    
    CMat = calculateSparseCoefficients(X,lambda,kernel);
elseif nargin<6
    CMat = calculateSparseCoefficients(X,lambda);
end

CKSym = BuildAdjacency(CMat,K);

[Grps , SingVals, LapKernel] = SpectralClustering(CKSym,n);
[Missrate, confusionMatrix,predicted] = Misclassification(Grps,groundTruth);

% display(confusionMatrix);
% display(min(Missrate));

end

