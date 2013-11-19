function [ alpha1 alphaStar ] = svrDual(kernel,y,C,epsilon )
%SVM Learn dual of soft-margin SVM with given data and paraemters and a
% given kernel
%
%   Input:
%
%       X               -       m x n data matrix or n x p cell array
%       K               -       kernel matrix of the function handle for
%           the kernel
%       y               -       y_i  for i=1,...,m
%       C               -       regularization parameter
%       epsilon         -       epsilon of the loss incensitive function
%
%   Output:
%
%       f               -       function handle to the classification
%           boundary
%
%  author: Ivan Bogun
%  date  : June 27, 2013

[n,~]=size(kernel);



cvx_begin quiet
variables alpha1(n) alphaStar(n)
        maximize (sum((alpha1-alphaStar).*y)-0.5* quad_form((alpha1-alphaStar)',kernel)-epsilon*sum(alpha1+alphaStar));
    subject to
        0<=alpha1;
        alpha1<=C;
        
        0<=alphaStar;
        alphaStar<=C;
        
        sum(alpha1-alphaStar)==0;
cvx_end