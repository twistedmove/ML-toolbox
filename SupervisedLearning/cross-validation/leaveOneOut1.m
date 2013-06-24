X=dataMatrix;

[n,m]=size(X);


classifierParameters=0.1;

y=-ones(m,1);

type=1;
X=X';

% ground truth
y(groundTruth==type)=1;


% create cross validation partition

results=zeros(54,6);
real=zeros(54,6);
for i=1:m
    fprintf('Current iteration %2g \n',i);
    
    
    for j=1:6
        type=j;
        y=-ones(m,1);
        y(groundTruth==type)=1;
        
        xTrain=X;
        yTrain=y;
        
        xTrain(i,:)=[];
        yTrain(i,:)=[];

        f=svm(xTrain,yTrain,classifierParameters);
        real(i,j)=svmPredict(f,X(i,:));
        %preds{i}=f;

        results(i,j)=f(X(i,:));

    end
end

R1=results-1; % 1
R2=results+1; % -1




%fprintf('%3g \n',sum(real==y)/m);