%% Project Overview Figures

% Producing figures for a project overview

clc, clear, close all

base = 'C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Project Overview Presentation\refs';

ff = '-dtiff'; %file format
p  = 1; %print? (aka save)

plot_where = [20,60];
plot_size  = [800,375];

%% Start the figure

figure('Position',[plot_where plot_size]), hold on
set(gca, 'FontSize', 16)
set(gcf,'defaultLineLineWidth',2)

%% SPDs

load sur_vrhel.mat
S_SRF = S_vrhel; clear S_vrhel
T_SRF = sur_vrhel(:,[1:44,65,69,81:154]);

xlim([390 730]); xticks(400:100:700);
ylim([0 1]); yticks([0,1]);
xlabel('Wavelength'); ylabel('Reflectance');

plot(SToWls(S_SRF),T_SRF(:,109)) %pear
if p
    print([base,'\','1'],ff)
end

%% all the SRFs
plot(SToWls(S_SRF),T_SRF)
if p
    print([base,'\','2'],ff)
end

%%

cla
[pc.coeff, pc.score, pc.latent, pc.tsquared, pc.explained, pc.mu] = pca(T_SRF');

ylim([-0.3 0.3]); yticks([-0.3,0,0.3]);
%legend({'PC1','PC2','PC3'},'Location','Southwest')
fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
    'k','LineStyle','none','FaceAlpha','0.03');

%plot(SToWls(S_SRF),pc.coeff(:,1:3)./max(pc.coeff(:,1:3)))

plot(SToWls(S_SRF),pc.coeff(:,1))
if p
    print([base,'\','3'],ff)
end
%%
plot(SToWls(S_SRF),pc.coeff(:,2))
if p
    print([base,'\','4'],ff)
end
%%
plot(SToWls(S_SRF),pc.coeff(:,3))
if p
    print([base,'\','5'],ff)
end

%%
cla
hold on
ylim([0 1]); yticks([0,1]);

recon(:,1) = pc.score(109,:)   * pc.coeff'        + pc.mu; % Reconstructed with all PCs
recon(:,2) = pc.score(109,1:3) * pc.coeff(:,1:3)' + pc.mu; % Reconstructed with 3 PCs

%plot(SToWls(S_SRF), T_SRF(:,109), 'b') % Original data
plot(SToWls(S_SRF), recon(:,2), 'r') % 3 PC
% legend({'Original Data','Reconstructed with 1:3'})
if p
    print([base,'\','6'],ff)
end

% figure, plot((pc.score*pc.coeff' + repmat(pc.mu,length(pc.score),1))') %all the data reconstructed

disp(pc.score(109,1:3))
