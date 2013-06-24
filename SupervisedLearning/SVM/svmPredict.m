function [ res ] = svmPredict( f,x )
%SVMPREDICT Summary of this function goes here
%   Detailed explanation goes here

y=f(x);

yRes=[y-1,y+1];
[~,indx]=min(abs(yRes));

if (indx==1)
    res=1;
else
    res=-1;
end

end

