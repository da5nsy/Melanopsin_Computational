clc, clear, close all

n = 7; %number of points per cluster
m = 4; %number of clusters

gt = repmat(randperm(m),n,1); %Ground Truth
gt = reshape(gt,1,[]);

kmi = gt(randperm(length(gt))); %fake K-Means-Index data (created by randomly shuffling gt)
%kmi = gt;

for i = 1:m
    mpc(i) = mean(mode(kmi(gt == i)) == kmi(gt == i)); %mode proportion
    % Issue: more than one group (defined by gt) could have the same mode. 
    % As m increases this should become less likely.
end

kmm = mean(mpc); %k-means-mark

% Where kmi is just a random permutation of gt the score is generally
% rather low (<0.5).
% Uncommenting line 8 (%kmi = gt;) shows what happens when the
% k-means-index performs perfectly: a score of 1.
