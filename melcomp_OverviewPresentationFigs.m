%% Project Overview Figures

% Producing figures for a project overview

clc, clear, close all

base = 'C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Project Overview Presentation';

ff = '-dtiff'; %file format
p  = 1; %print? (aka save?), set to 1 to begin saving

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

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% all the SPDs
plot(SToWls(S_SPD),T_SPD)
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%%

cla
[pc_p.coeff, pc_p.score, pc_p.latent, pc_p.tsquared, pc_p.explained, pc_p.mu] = pca(T_SPD');

ylim([-0.3 0.3]); yticks([-0.3,0,0.3]);
%legend({'PC1','PC2','PC3'},'Location','Southwest')
fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
    'k','LineStyle','none','FaceAlpha','0.03');

%plot(SToWls(S_SPD),pc.coeff(:,1:3)./max(pc.coeff(:,1:3)))

plot(SToWls(S_SPD),pc_p.coeff(:,1))
if p
    print([base,'\','3'],ff)
end
%%
plot(SToWls(S_SPD),pc_p.coeff(:,2))
if p
    print([base,'\','4'],ff)
end
%%
plot(SToWls(S_SPD),pc_p.coeff(:,3))
if p
    print([base,'\','5'],ff)
end

%%
cla
hold on
ylim([0 2]); yticks([0,1,2]);

recon_p(:,1) = pc_p.score(103,:)   * pc_p.coeff'        + pc_p.mu; % Reconstructed with all PCs
recon_p(:,2) = pc_p.score(103,1:3) * pc_p.coeff(:,1:3)' + pc_p.mu; % Reconstructed with 3 PCs

%plot(SToWls(S_SPD), T_SPD(:,103), 'b') % Original data
plot(SToWls(S_SPD), recon_p(:,2), 'r') % 3 PC
% legend({'Original Data','Reconstructed with 1:3'})
if p
    print([base,'\','6'],ff)
end

% figure, plot((pc_p.score*pc_p.coeff' + repmat(pc_p.mu,length(pc_p.score),1))') %all the data reconstructed

%% ----------------------%

%% SRFs

load sur_vrhel.mat
S_SRF = S_vrhel; clear S_vrhel
T_SRF = sur_vrhel(:,[1:44,65,69,81:154]);

cla
xlim([390 730]); xticks(400:100:700);
ylim([0 1]); yticks([0,1]);
xlabel('Wavelength'); ylabel('Reflectance');

plot(SToWls(S_SRF),T_SRF(:,109)) %pear
if p
    print([base,'\','7'],ff)
end

%% all the SRFs
plot(SToWls(S_SRF),T_SRF)
if p
    print([base,'\','8'],ff)
end

%%

cla
[pc_r.coeff, pc_r.score, pc_r.latent, pc_r.tsquared, pc_r.explained, pc_r.mu] = pca(T_SRF');

ylim([-0.3 0.3]); yticks([-0.3,0,0.3]);
%legend({'PC1','PC2','PC3'},'Location','Southwest')
fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
    'k','LineStyle','none','FaceAlpha','0.03');

%plot(SToWls(S_SRF),pc.coeff(:,1:3)./max(pc.coeff(:,1:3)))

plot(SToWls(S_SRF),pc_r.coeff(:,1))
if p
    print([base,'\','9'],ff)
end
%%
plot(SToWls(S_SRF),pc_r.coeff(:,2))
if p
    print([base,'\','10'],ff)
end
%%
plot(SToWls(S_SRF),pc_r.coeff(:,3))
if p
    print([base,'\','11'],ff)
end

%%
cla
hold on
ylim([0 1]); yticks([0,1]);

recon_r(:,1) = pc_r.score(109,:)   * pc_r.coeff'        + pc_r.mu; % Reconstructed with all PCs
recon_r(:,2) = pc_r.score(109,1:3) * pc_r.coeff(:,1:3)' + pc_r.mu; % Reconstructed with 3 PCs

%plot(SToWls(S_SRF), T_SRF(:,109), 'b') % Original data
plot(SToWls(S_SRF), recon_r(:,2), 'r') % 3 PC
% legend({'Original Data','Reconstructed with 1:3'})
if p
    print([base,'\','12'],ff)
end

% figure, plot((pc_r.score*pc_r.coeff' + repmat(pc_r.mu,length(pc_r.score),1))') %all the data reconstructed

% disp(pc_r.score(109,1:3)) %values for demo ref

%% Plot spectral sensitivities

figure(1)
cla

load T_cones_sp T_cones_sp S_cones_sp
load T_rods T_rods S_rods
load T_melanopsin T_melanopsin S_melanopsin

plot(SToWls(S_cones_sp),T_cones_sp)
plot(SToWls(S_rods),T_rods)
plot(SToWls(S_melanopsin),T_melanopsin)

ylabel('Sensitivity')

if p
    print([base,'\','13'],ff)
end

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

origin.score = median(pc_p.score(:,1:3));
origin.spd   = origin.score * pc_p.coeff(:,1:3)' + pc_p.mu; %plot(origin.spd)
origin.XYZ   = T_xyz1931 * origin.spd';
origin.xy    = [origin.XYZ(1)./sum(origin.XYZ);origin.XYZ(2)./sum(origin.XYZ)];
%scatter(origin.xy(1),origin.xy(2),'k*')

for i=1:3
    n.incpc(:,i) = [linspace(min(pc_p.score(:,i)),median(pc_p.score(:,i)),20),...
        linspace(median(pc_p.score(:,i)),max(pc_p.score(:,i)),20)];
    % - incremental pc weights
    % - Should find a way to remove middle repeating value
    % - 5th/95th percentile might be preferable to min/max
end
%figure, plot(n.incpc)

n.score = repmat(origin.score,length(n.incpc),1,3);
%n.score = zeros(40,3,3);
for i=1:3
    n.score(:,i,i) = n.incpc(:,i);
end

n.spd = zeros(length(n.incpc),3,size(pc_p.score,2));

for i=1:length(n.incpc)
    for j=1:3
        n.spd(i,j,:) = n.score(i,:,j) * pc_p.coeff(:,1:3)' + pc_p.mu;
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
    print([base,'\','14'],ff)
end
scatter3(n.xy(:,2,1),n.xy(:,2,2),n.XYZ(:,2,2),'b*')
if p
    print([base,'\','15'],ff)
end
scatter3(n.xy(:,3,1),n.xy(:,3,2),n.XYZ(:,3,2),'r*')
if p
    print([base,'\','16'],ff)
end

%%
figure(1)
cols = lines(3);

cla
ylim([-0.3 0.3]); yticks([-0.3,0,0.3]);
fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
    'k','LineStyle','none','FaceAlpha','0.03');

plot(SToWls(S_SPD),pc_p.coeff(:,2),'Color',cols(2,:))
if p
    print([base,'\','17'],ff)
end

for i = -2:0.5:2
    plot(SToWls(S_SPD),pc_p.coeff(:,2)*i,'Color',cols(2,:)*((i+2)/4))
end
if p
    print([base,'\','18'],ff)
end

cla
plot(SToWls(S_SPD),pc_p.coeff(:,2)*-1,'Color',cols(2,:)*((-1+2)/4))
if p
    print([base,'\','19'],ff)
end

%%

