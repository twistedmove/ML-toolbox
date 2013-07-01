X=trajectoriesArray(1,:);
n=length(X);
y=-ones(1,54);

% current class learning
k=1;
C=1;
kernel=getKernel('gaussianAlignment',10);
y(groundTruth==k)=1;
f=trainSVM( X,y',kernel,C );