

y=groundTruth;

n=length(groundTruth);

y=-ones(n,1);
k=1;
y(groundTruth==k)=1;


kernel{1}=getKernel('gaussianAlignment',100);
kernel{2}=getKernel('gaussianAlignment',10);

K=calculateKernel(trajectoriesArray,kernel);

weights1=getWeightsMKL(K,y,'QuiLane2009');
weights2=getWeightsMKL(K,y,'Lanckriet2004a');
weights3=getWeightsMKL(K,y,'Cortes2010a');

