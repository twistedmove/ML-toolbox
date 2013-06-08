function [ kernel ] = getKernel( kernelType, params )
%GETKERNEL The function will return a kernel function handle
%   The function will return the kernel specified by the string 'kernelType' 
%       Parameters of the kernel are given by the verctor params.
%       Input for the kernels have to be a column vectors.
%
%   Input: 
%       kernelType          -       name of the kernel 
%       params              -       parameters of the kernel
%
%   Output:
%       kernel              -       function handle of the kernel
%
%   author: Ivan Bogun
%   date  : June 8, 2013
%

if nargin<2
    params=1;
end


% there is no such function ( only for convexity of certain kernels in CVX)
if (exist('pow_pos')==0)
    pow_pos=@(x,p) max(x^p,0);
end


if (strcmp(kernelType,'linear'))
    kernel=@(x,y) x'*y;
    
elseif(strcmp(kernelType,'polynomial'))
    kernel=@(x,y) (x'*y +1)^params(1);
    
elseif(strcmp(kernelType,'gaussian'))
    kernel=@(x,y) exp( -(pow_pos(norm(x-y),2))/(2*params(1)^2));
    
elseif(strcmp(kernelType,'laplacian'))
    kernel=@(x,y) exp( -((norm(x-y)))/(params(1)));
    
elseif(strcmp(kernelType,'laplacian'))
    kernel=@(x,y) exp( -((norm(x-y)))/(2*params(1)));
    
elseif(strcmp(kernelType,'rationalQuadratic'))
    kernel=@(x,y) 1 - (pow_pos(norm(x-y),2)) /(pow_pos(norm(x-y),2)+params(1));  
        
elseif(strcmp(kernelType,'multiquadratic'))
    kernel=@(x,y) sqrt( pow_pos(norm(x-y),2)+params(1)^2);  
            
elseif(strcmp(kernelType,'inverseMultiQuadratic'))
    kernel=@(x,y) 1/sqrt( pow_pos(norm(x-y),2)+params(1)^2); 
            
elseif(strcmp(kernelType,'wave'))
    kernel=@(x,y) (params(1)/norm(x-y)) * sin( norm(x-y)/params(1)); 
            
elseif(strcmp(kernelType,'cauchy'))
    kernel=@(x,y) 1/(1+ pow_pos(norm(x-y),2)/params(1)); 
            
elseif(strcmp(kernelType,'chi-square'))
    kernel=@(x,y) sum((x-y).^2./(x+y)/2); 
                
elseif(strcmp(kernelType,'generalized_T-student'))
    kernel=@(x,y) 1/(1+pow_p(norm(x-y),params(1))); 
    
end


end

