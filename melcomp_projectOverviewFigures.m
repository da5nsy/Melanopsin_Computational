%% Project Overview Figures

% Producing figures for a project overview

clc, clear, close all

base = 'C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Project Overview Presentation';

ff = '-dtiff'; %file format
p  = 0; %print? (aka save)

plot_where = [20,60];
plot_size  = [800,375];

%% Start the figure

figure('Position',[plot_where plot_size]), hold on
set(gca, 'FontSize', 16)
set(gcf,'defaultLineLineWidth',2)

%% SPDs

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)

T_SPD = final(17:97,:); clear final
S_SPD = [380,5,81];

xlim([380 780]); xticks(380:100:780);
ylim([0 2]); yticks([0,1,2]);
xlabel('Wavelength'); ylabel('Power');

plot(SToWls(S_SPD),T_SPD(:,103))
if p
    print([base,'\','1'],ff)
end

%% all the SPDs
plot(SToWls(S_SPD),T_SPD)
if p
    print([base,'\','2'],ff)
end

%%

cla
[pc.coeff, pc.score, pc.latent, pc.tsquared, pc.explained, pc.mu] = pca(T_SPD');

ylim([-0.3 0.3]); yticks([-0.3,0,0.3]);
%legend({'PC1','PC2','PC3'},'Location','Southwest')
fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
    'k','LineStyle','none','FaceAlpha','0.03');

%plot(SToWls(S_SPD),pc.coeff(:,1:3)./max(pc.coeff(:,1:3)))

plot(SToWls(S_SPD),pc.coeff(:,1))
if p
    print([base,'\','3'],ff)
end
%%
plot(SToWls(S_SPD),pc.coeff(:,2))
if p
    print([base,'\','4'],ff)
end
%%
plot(SToWls(S_SPD),pc.coeff(:,3))
if p
    print([base,'\','5'],ff)
end

%%
cla
hold on
ylim([0 2]); yticks([0,1,2]);

recon(:,1) = pc.score(103,:)   * pc.coeff'        + pc.mu; % Reconstructed with all PCs
recon(:,2) = pc.score(103,1:3) * pc.coeff(:,1:3)' + pc.mu; % Reconstructed with 3 PCs

%plot(SToWls(S_SPD), T_SPD(:,103), 'b') % Original data
plot(SToWls(S_SPD), recon(:,2), 'r') % 3 PC
% legend({'Original Data','Reconstructed with 1:3'})
if p
    print([base,'\','6'],ff)
end

% figure, plot((pc.score*pc.coeff' + repmat(pc.mu,2600,1))') %all the data reconstructed



