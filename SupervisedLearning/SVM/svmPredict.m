function [ label value] = svmPredict( f,x )
%SVMPREDICT Summary of this function goes here
%   Detailed explanation goes here

value=f(x);

yRes=[value-1,value+1];
[~,indx]=min(abs(yRes));

if (indx==1)
    label=1;
else
    label=-1;
end

end

