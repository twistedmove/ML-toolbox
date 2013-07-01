function [ weights ] = getWeightsMKL( K,y,method )
%GETWEIGHTSMKL This function will learn weights for the MKL problem.
%
%   Input:
%       K           -       NxNxM matrix which represents M kernels
%       y           -       Mx1 label matrix with y\in {1,-1}
%
%   Output:
%       weights     -       Mx1 matrix with kernel weights
%
% Methods implemented: 'QuiLane2009','Lanckriet2004a', 'Cortes2010a'
% from the paper:
%   "Multiple kernel Learning Algorithms" by Mehmet Gonen, Ethem Alpaydin
%
%   author: Ivan Bogun
%   date  : July 1, 2013


idealKernel=y*y';

m=size(K,3);
weights=zeros(m,1);

switch method
    case 'QuiLane2009'
        
        for i=1:m
            weights(i,1)=kernelAlignment(K(:,:,i),idealKernel);
            
        end
        weights=weights/sum(weights);
        
    case 'Lanckriet2004a'
                
        M=zeros(m,m);
        a=zeros(m,1);
        
        for i=1:m
            
            varK=(K(:,:,i));
            a(i,1)=kernelAlignment(varK,idealKernel);
            
            for j=1:m
                
                M(i,j)=kernelAlignment(varK,(K(:,:,j)));
                
            end
        end
        
        
        cvx_begin quiet
        
        variable nu(m)
        
        maximize (nu'*a)
        subject to
        quad_form(nu,M)<=1
        nu>=0
        
        cvx_end
        
        weights=nu;
        
    case 'Cortes2010a'
        
        M=zeros(m,m);
        a=zeros(m,1);
        
        for i=1:m
            
            varK=centerKernelAlignment(K(:,:,i));
            a(i,1)=kernelAlignment(varK,idealKernel);
            
            for j=1:m
                
                M(i,j)=kernelAlignment(varK,...
                    centerKernelAlignment(K(:,:,j)));
                
                
                
            end
        end
        
        cvx_begin quiet
        variables v(m)
        minimize (quad_form(v,M)-2*v'*a)
        subject to
        v>=0;
        cvx_end
        
        weights=v/norm(v);
    otherwise
        display('Not supported');
        
        
end

end


function res=kernelAlignment(K1,K2)

inner=@(x,y) x(:)'*y(:);

res=inner(K1,K2)/(sqrt(inner(K1,K1)*inner(K2,K2)));

end

function Kc=centerKernelAlignment(K)

n=size(K,1);
one=ones(n,1);

Kc=(eye(n)-(one*one')/n)*K*(eye(n)-(one*one')/n);

end
