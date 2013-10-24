function [ CMat ] = calculateSparseCoefficients( Xp,lambda,kernel)
%CALCULATESPARSECOEFFICIENTS Calculate sparse regressors as a part of SSC
%algorithm
%   Convex optimization problem will be solved for sparse regressors of the
%   SSC algorithm
%
%   Input:
%       Xp          -       DxN data matrix of N data points
%       lambda      -       regularization parameter
%       kernel      -       kernel to use in the kernel version of SSC, can
%       be a struct so that kernel(i).kernel - function handle,
%       kernel(i).weights - weight of the kernel i.
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


if (iscell(Xp)==1)
    fullK=kernel;
    %     if (size(Xp,1)==1)
    %
    %         for s=1:(N)
    %             for j=1:(N)
    %                 if (j>s)
    %                     fullK(s,j)=kernel(Xp{1,s},Xp{1,j});
    %                 elseif (j==s)
    %                     fullK(j,s)=1;
    %                 else
    %                     fullK(s,j)=fullK(j,s);
    %                 end
    %             end
    %         end
    %
    %     else
    %
    %         for k=1:size(Xp,1)
    %             fullK1=zeros(N,N);
    %
    %             for s=1:(N)
    %                 for j=1:(N)
    %                     if (j>=s)
    %                         fullK1(s,j)=kernel(k).kernels(Xp{2,s},Xp{2,j});
    %                     else
    %                         fullK1(s,j)=fullK1(j,s);
    %                     end
    %                 end
    %             end
    %             fullK=fullK+kernel(k).weights*fullK1;
    %         end
    %
    %     end
    %     minEigenvalue=min(eig(fullK));
    %     if (minEigenvalue<0)
    %         fullK=fullK-minEigenvalue*eye(N);
    %     end
end


for i = 1:N
    
    
    %fprintf('%3g data record out of %3g \n',i,N);
    
    if (iscell(Xp)==1)
        y=Xp{1,i};
    else
        y = Xp(:,i);
    end
    
    
    if i == 1
        Y = Xp(:,i+1:end);
    elseif ( (i > 1) && (i < N) )
        Y = [Xp(:,1:i-1) Xp(:,i+1:N)];
    else
        Y = Xp(:,1:N-1);
    end
    
    
    if (kernelMode==0)
        
        % original SSC algorithm
        if (iscell(Xp)==1)
            error('Linear SSC is not implemented for cell arrays, use kernel version instead');
        else
            cvx_begin;
            
            cvx_quiet true;
            cvx_precision high
            
            variable c(N-1,1);
            
            minimize( norm(c,1) + lambda * norm(Y * c  - y) );
            
            subject to
            
            cvx_end;
        end
    else
        
        % kernel variant of the SSC algorithm
        
        
        if (iscell(Xp)==1)
            
            K_yy=fullK(i,i);
            
            Kv_y=fullK(:,i);
            Kv_y(i)=[];
            
            K=fullK;
            K(i,:)=[];
            K(:,i)=[];
            
        else
            
            for s=1:(N-1)
                Kv_y(s)=kernel(Y(:,s),y);
                for j=1:(N-1)
                    K(s,j)=kernel(Y(:,s),Y(:,j));
                end
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

