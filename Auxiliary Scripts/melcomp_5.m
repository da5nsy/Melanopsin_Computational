clc, clear, close all

%% load granada

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
T_SPD=final; clear final
S_SPD=[300,5,161];

%calculate linear signals and mb chromaticities for daylights
%   plots

%% calculate pca of daylight
[pc_p.coeff, pc_p.score, pc_p.latent, pc_p.tsquared, pc_p.explained, pc_p.mu] = pca(T_SPD');

%create set of illuminants with median pc1 + pc2
recon_p = repmat(median(pc_p.score(:,1)),2600,1) * pc_p.coeff(:,1)' + pc_p.score(:,2) * pc_p.coeff(:,2)';

recon_p2 = repmat(median(pc_p.score(:,1)),2600,1) * pc_p.coeff(:,1)' + pc_p.score(:,2:end) * pc_p.coeff(:,2:end)';

recon_p3 = repmat(median(pc_p.score(:,1)),2600,1) * pc_p.coeff(:,1)' + pc_p.score(:,3) * pc_p.coeff(:,3)';

recon_p3 = pc_p.score(:,1) * pc_p.coeff(:,1)' + pc_p.score(:,2:end) * pc_p.coeff(:,2:end)' + pc_p.mu;
%create set of illuminants with median pc1 + linspace(min(pc2,max(pc2),100(?))

%create signals for made-up cone types
%   add noise
%   comp the correlation is between each signal and the value of pc2

figure,
plot(SToWls(S_SPD),recon_p')
figure,
plot(SToWls(S_SPD),recon_p2')
figure,
plot(SToWls(S_SPD),recon_p3')

figure,
plot(SToWls(S_SPD),std(recon_p))

%%

