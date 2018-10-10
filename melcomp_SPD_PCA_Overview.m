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
[pc.coeff, pc.score, pc.latent, pc.tsquared, pc.explained, pc.mu] = pca(T_SPD','Centered',false);

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

% figure, plot((pc.score*pc.coeff' + repmat(pc.mu,length(pc.score),1))') %all the data reconstructed

%% Show colorimetric effect of principal components

%load CMF
load T_xyz1931.mat

%start figure
figure('Position',[plot_where 800 800]), hold on
set(gca, 'FontSize', 16)
set(gcf,'defaultLineLineWidth',2)
axis square
xlim([0 1]); xticks([0,1]);
ylim([0 1]); yticks([0,1]);
xlabel('x')
ylabel('y')

%plot spectral locus
plot(T_xyz1931(1,:)./sum(T_xyz1931),T_xyz1931(2,:)./sum(T_xyz1931),...
    'ko-','MarkerFaceColor','k','MarkerSize',2);

origin.score = median(pc.score(:,1:3));
origin.spd   = origin.score * pc.coeff(:,1:3)' + pc.mu; %plot(origin.spd)
origin.XYZ   = T_xyz1931 * origin.spd';
origin.xy    = [origin.XYZ(1)./sum(origin.XYZ);origin.XYZ(2)./sum(origin.XYZ)];
%scatter(origin.xy(1),origin.xy(2),'k*')

for i=1:3
    n.incpc(:,i) = [linspace(min(pc.score(:,i)),median(pc.score(:,i)),20),...
        linspace(median(pc.score(:,i)),max(pc.score(:,i)),20)];
    % - incremental pc weights
    % - Should find a way to remove middle repeating value
    % - 5th/95th percentile might be preferable to min/max
end
%figure, plot(n.incpc)

%n.score = repmat(origin.score,length(n.incpc),1,3);
n.score = zeros(40,3,3);
for i=1:3
    n.score(:,i,i) = n.incpc(:,i);
end

n.spd = zeros(length(n.incpc),3,size(pc.score,2));

for i=1:length(n.incpc)
    for j=1:3
        n.spd(i,j,:) = n.score(i,:,j) * pc.coeff(:,1:3)' + pc.mu;
    end
end

% figure, hold on, axis tight
% plot(SToWls(S_SPD),squeeze(n.spd(:,1,:))'/max(max(squeeze(n.spd(:,1,:)))),'k')
% plot(SToWls(S_SPD),squeeze(n.spd(:,2,:))'/max(max(squeeze(n.spd(:,2,:)))),'b')
% plot(SToWls(S_SPD),squeeze(n.spd(:,3,:))'/max(max(squeeze(n.spd(:,3,:)))),'r')

for i=1:length(n.incpc)
    for j=1:3
        n.XYZ(i,j,:) = T_xyz1931 * squeeze(n.spd(i,j,:));
        n.xy(i,j,:)  = [n.XYZ(i,j,1)./sum(n.XYZ(i,j,:));n.XYZ(i,j,2)./sum(n.XYZ(i,j,:))];
    end
end

scatter3(n.xy(:,1,1),n.xy(:,1,2),n.XYZ(:,1,2),'k*')
zlabel('Y')
if p
    print([base,'\','7'],ff)
end
scatter3(n.xy(:,2,1),n.xy(:,2,2),n.XYZ(:,2,2),'b*')
if p
    print([base,'\','8'],ff)
end
scatter3(n.xy(:,3,1),n.xy(:,3,2),n.XYZ(:,3,2),'r*')
if p
    print([base,'\','9'],ff)
end

