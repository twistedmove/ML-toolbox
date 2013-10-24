

% number of kernels
n=2;


kernels{1}=getKernel('gaussianAlignment',1);
kernels{2}=getKernel('gaussianAlignment',1);

res=MKSSC(trajectoriesArray,6,0,10,kernels);


% CKSym = BuildAdjacency(CMat,0);
% 
% [Grps , ~, ~] = SpectralClustering(CKSym,6);
% load('groundTruth');
% [Missrate, confusionMatrix,predicted] = Misclassification(Grps,groundTruth);
% display(min(Missrate));