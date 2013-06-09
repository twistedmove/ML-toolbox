function [ CMat ] = calculateSparseCoefficients( Xp,lambda,kernel)
%CALCULATESPARSECOEFFICIENTS Calculate sparse regressors as a part of SSC
%algorithm
%   Convex optimization problem will be solved for sparse regressors of the
%   SSC algorithm
%
%   Input:
%       Xp          -       DxN data matrix of N data points
%       lambda      -       regularization parameter
%       kernel      -       kernel to use in the kernel version of SSC
%
%   Output:
%       CMat        -       regressor coefficients
%
%   author: Ivan Bogun
%   date  : June 8, 2013
%
%   Credit: the code was adapted from the one by Ehsan Elhamifar


N = size(Xp,2);

if (nargin < 2)
    lambda = 0.001;
end

if nargin==3
    kernelMode=1;
    
    K=zeros(N-1,N-1);
    Kv_y=zeros(N-1,1);
else
    kernelMode=0;
end




for i = 1:N
    
    fprintf(1,'%3g data record out of %3g \n',i,N);
    
    y = Xp(:,i);
    if i == 1
        Y = Xp(:,i+1:end);
    elseif ( (i > 1) && (i < N) )
        Y = [Xp(:,1:i-1) Xp(:,i+1:N)];
    else
        Y = Xp(:,1:N-1);
    end
    
    
    if (kernelMode==0)
        
        % original SSC algorithm
        cvx_begin;
        
        cvx_quiet true;
        cvx_precision high
        
        variable c(N-1,1);
        
        minimize( norm(c,1) + lambda * norm(Y * c  - y) );
        
        subject to
        
        cvx_end;
        
    else
        
        % kernel variant of the SSC algorithm
        
        K_yy=kernel(y,y);
        
        for s=1:(N-1)
            Kv_y(s)=kernel(Y(:,s),y);
            for j=1:(N-1)
                K(s,j)=kernel(Y(:,s),Y(:,j));
            end
        end
        
        %version 1
        cvx_begin;
        
        cvx_quiet true;
        cvx_precision high
        
        variable c(N-1,1);
        
        minimize( norm(c,1)+lambda*(K_yy-2*c'*(Kv_y)+quad_form(c,K)));
        
        subject to
        
        cvx_end;
        
        %   version 2
        %     cvx_begin;
        %
        %     cvx_quiet true;
        %     cvx_precision high
        %     variable c(N-1,1);
        %     minimize( norm(c,1));
        %     subject to
        %     (K_yy-2*c'*(Kv_y)+quad_form(c,K))<=lambda;
        %     cvx_end;
        
    end
    
    
    % place 0's in the diagonals of the coefficient matrix
    if i == 1
        CMat(1,1) = 0;
        CMat(2:N,1) = c(1:N-1);
    elseif ( (i > 1) && (i < N) )
        CMat(1:i-1,i) = c(1:i-1);
        CMat(i,i) = 0;
        CMat(i+1:N,i) = c(i:N-1);
    else
        CMat(1:N-1,N) = c(1:N-1);
        CMat(N,N) = 0;
    end
    
    
end

