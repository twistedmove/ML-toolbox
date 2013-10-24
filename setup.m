host=pwd;

addpath_recurse(strcat(host,'/Other')); 
addpath_recurse(strcat(host,'/Probabilistic')); 
addpath_recurse(strcat(host,'/SupervisedLearning')); 
addpath_recurse(strcat(host,'/UnsupervisedLearning')); 

addpath_recurse('Other'); 
addpath_recurse('Probabilistic'); 
addpath_recurse('SupervisedLearning'); 
addpath_recurse('UnsupervisedLearning');
clearvars host;