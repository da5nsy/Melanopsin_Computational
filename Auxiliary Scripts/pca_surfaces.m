clear, clc, close all

load sur_nickerson.mat

figure, 
plot(SToWls(S_nickerson),sur_nickerson)


%%

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(sur_nickerson');

figure,
plot(COEFF(:,1:4))

legend

%%
load B_cohen.mat

figure, 
plot(B_cohen)

legend

%%

% 'helix' shape noted by Cohen (1964)
figure,
plot3(COEFF(:,1),COEFF(:,2),COEFF(:,3),'.')



