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

n=8;

figure, 
plot(SToWls(S_vrhel),coeff(:,1:n).*explained(1:n)','DisplayName',['PC',num2str(n)])
legend('Location','best')

xlabel('Wavelength (nm)')
ylabel('coeff * explained')