function [ f ] = svm(X,y,C)
%SVM Learn soft-margin SVM with given data and paraemters
%   
%   Input:
%
%       X               -       [m,n] data matrix with m data points each
%           with dimension n
%       y               -       y_i in {-1,1} for i=1,...,m
%       C               -       regularization parameter
%
%   Output:
%
%       f               -       function handle to the classification
%           boundary
%
%  author: Ivan Bogun
%  date  : June 19, 2013

% TO DO: solve kernel dual instead of linear SVM

[m,n]=size(X);

cvx_begin quiet
    variables w(n) b err(m)
    minimize 1/2*sum(w.*w) + C*sum(err)
    subject to
        y.*(X*w + b) >= 1 - err;
        err >= 0;
cvx_end

f=@(x) x*w+b;

end

