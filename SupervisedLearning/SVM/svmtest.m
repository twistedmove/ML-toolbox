

n=100;
m=50;

y=ones(m,1);
y(25:m)=-1;

X=rand(m,n);
C=100;
% % m data points
% % n size of the each vector
% cvx_begin quiet
%     variables w(n) b err(m)
%     minimize 1/2*sum(w.*w) + C*sum(err)
%     subject to
%         y.*(X*w + b) >= 1 - err;
%         err >= 0;
% cvx_end

f=svm(X,y,C);