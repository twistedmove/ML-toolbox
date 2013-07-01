function [ alpha ] = svmDual(kernel,y,C )
%SVM Learn dual of soft-margin SVM with given data and paraemters and a
% given kernel
%
%   Input:
%
%       X               -       m x n data matrix or n x p cell array
%       K               -       kernel matrix of the function handle for
%           the kernel
%       y               -       y_i in {-1,1} for i=1,...,m
%       C               -       regularization parameter
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
variables alpha(n)
        maximize (sum(alpha) -0.5* quad_form(alpha.*y,kernel))
    subject to
        0<=alpha;
        alpha<=C;
        alpha'*y==0;
cvx_end

% y_neg=find(y==-1);
% y_pos=find(y==1);
% 
% w=X*(alpha.*y);
% 
% yMax=zeros(length(y_neg),1);
% yMin=zeros(length(y_pos),1);
% 
% % for i=1:length(y_neg)
% %     yMax(i)=w'*
% % end
% b=0;
% f=@(x) x*w+b;




end

