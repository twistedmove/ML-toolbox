%   pLSA test cross Validation
clc; close all; clear all;

load('pLSA_Testing.mat');

videos = length(groundTruth);

predictedLabels = zeros(videos, 1);

for video = 1 : videos
    labels = groundTruth;
    labels(video) = 0;
    [Pw_z,Pz_d,Pz,Li] = pLSA_EMmodified(dataMatrix, numberOfClasses, Par,labels);
    [~, predictedLabels(video)] = max(Pz_d(video, :));
end
labels = groundTruth;
sum(predictedLabels==groundTruth)