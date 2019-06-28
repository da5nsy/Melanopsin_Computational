clear, clc, close all

%%

load sur_vrhel.mat

[coeff, score, latent, tsquared, explained, mu] = pca(sur_vrhel');

%%

figure,
scatter(1:10,explained(1:10),'k','filled')
ylim([0 100])
grid on

%%

figure, hold on
plot(SToWls(S_vrhel),coeff(:,1))
plot(SToWls(S_vrhel),coeff(:,2))
plot(SToWls(S_vrhel),coeff(:,3))
axis tight

legend('Location','best')