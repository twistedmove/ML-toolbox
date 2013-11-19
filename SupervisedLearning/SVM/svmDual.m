function [ alpha1 ] = svmDual(kernel,y,C,classPrior)
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

if nargin>3
   % we know some class prior ( data is not distributed evenly)
   
   pos=1-(sum(y==1))/length(y);
   neg=1-(sum(y==-1))/length(y);
   
   prior=zeros(length(y),1);
   
   posIdx=find(y==1,length(y),'first');
   negIdx=find(y==-1,length(y),'first');
   
   prior(posIdx)=pos;
   prior(negIdx)=neg;
else
    
    prior=1;
    
end

cvx_begin quiet
variables alpha1(n)
        maximize (sum(alpha1) -0.5* quad_form(alpha1.*y,kernel))
    subject to
        0<=alpha1;
        alpha1<=C*prior;
        alpha1'*y==0;
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

