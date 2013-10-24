function [ bestAlignment ] = alignment( A,B,loss )
%ALIGNMENT Align sequences A and B according to the loss function
%
%   Input: 
%       A                 -       k x n sequence
%       B                 -       k x m sequence
%       loss              -       loss function
%
%   Ouput:
%       bestAlignment     -   score of the best alignment of the sequences
%           according to the given loss function
%
%
%   author: Ivan Bogun
%   date  : June 22, 2013

if nargin<3
    loss=@(a,b) -norm(a-b);
end

[n,k]=size(A);
m=size(B,1);

M=zeros(n+1,m+1);

M(1,1)=0;
null=zeros(1,k);
for i=2:n+1
    
    M(i,1)=M(i-1,1)+loss(A(i-1,:),null);
  
end

for j=2:m+1
    M(1,j)=M(1,j-1)+loss(null,B(j-1,:));
end

for i=2:n+1
    for j=2:m+1
        
        M(i,j)=max([M(i-1,j-1)+loss(A(i-1,:),B(j-1,:)),...
            M(i-1,j)+loss(A(i-1,:),null),...
            M(i,j-1)+loss(null,B(j-1,:))]);
    end
end
bestAlignment=M(end,end);
end

