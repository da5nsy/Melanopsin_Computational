% PCA of vrhel spectral reflectances

% Reproducing figure 20.5 of:
% Chapter 20: Physics-based approaches to modeling surface color perception. Laurence T. Maloney
% in:
% K. R. Gegenfurtner and L. T. Sharpe, Color Vision: From Genes to Perception. Cambridge University Press, 2001.

clear, clc, close all

load sur_vrhel.mat

p = pca(sur_vrhel');

figure, hold on
for i=1:8
    subplot(2,4,i)
    plot(SToWls(S_vrhel),p(:,i))
    axis tight
    %ylim([-1 1])
end

% Seems pretty similar
% Couple of inversions
% Differences at high wavelengths for 1 and 3

%% Nat only

sur_vrhel_n = sur_vrhel(:,[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141]);
plot(sur_vrhel_n)

p_n = pca(sur_vrhel_n');

figure, hold on
for i=1:8
    subplot(2,4,i)
    plot(SToWls(S_vrhel),p_n(:,i))
    axis tight
    %ylim([-1 1])
end

% Less neat but what did I expect with that few refs



