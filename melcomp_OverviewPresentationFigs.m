%% Project Overview Figures

% Producing figures for a project overview

% Now merged into melcomp_3.m

clc, clear, close all

base = 'C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Project Overview Presentation';

ff = '-dtiff'; %file format
p  = 0; %print? (aka save?), set to 1 to begin saving

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
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
%%
plot(SToWls(S_SPD),pc_p.coeff(:,2))
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
%%
plot(SToWls(S_SPD),pc_p.coeff(:,3))
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%%
cla
hold on
ylim([0 2]); yticks([0,1,2]);

recon_p(:,1) = pc_p.score(103,:)   * pc_p.coeff'        + pc_p.mu; % Reconstructed with all PCs
recon_p(:,2) = pc_p.score(103,1:3) * pc_p.coeff(:,1:3)' + pc_p.mu; % Reconstructed with 3 PCs

%plot(SToWls(S_SPD), T_SPD(:,103), 'b') % Original data
plot(SToWls(S_SPD), recon_p(:,2), 'r') % 3 PC
% legend({'Original Data','Reconstructed with 1:3'})
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

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
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% all the SRFs
plot(SToWls(S_SRF),T_SRF)
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%%

cla
[pc_r.coeff, pc_r.score, pc_r.latent, pc_r.tsquared, pc_r.explained, pc_r.mu] = pca(T_SRF');

ylim([-0.3 0.3]); yticks([-0.3,0,0.3]);
%legend({'PC1','PC2','PC3'},'Location','Southwest')
fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
    'k','LineStyle','none','FaceAlpha','0.03');

%plot(SToWls(S_SRF),pc.coeff(:,1:3)./max(pc.coeff(:,1:3)))

plot(SToWls(S_SRF),pc_r.coeff(:,1))
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
%%
plot(SToWls(S_SRF),pc_r.coeff(:,2))
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
%%
plot(SToWls(S_SRF),pc_r.coeff(:,3))
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%%
cla
hold on
ylim([0 1]); yticks([0,1]);

recon_r(:,1) = pc_r.score(109,:)   * pc_r.coeff'        + pc_r.mu; % Reconstructed with all PCs
recon_r(:,2) = pc_r.score(109,1:3) * pc_r.coeff(:,1:3)' + pc_r.mu; % Reconstructed with 3 PCs

%plot(SToWls(S_SRF), T_SRF(:,109), 'b') % Original data
plot(SToWls(S_SRF), recon_r(:,2), 'r') % 3 PC
% legend({'Original Data','Reconstructed with 1:3'})
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

% figure, plot((pc_r.score*pc_r.coeff' + repmat(pc_r.mu,length(pc_r.score),1))') %all the data reconstructed

% disp(pc_r.score(109,1:3)) %values for demo ref

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
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

scatter3(n.xy(:,2,1),n.xy(:,2,2),n.XYZ(:,2,2),'b*') %could be nice to use 'real' colours here(?)
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

scatter3(n.xy(:,3,1),n.xy(:,3,2),n.XYZ(:,3,2),'r*')
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% Plot spectral sensitivities

figure(1)
cla

ylim([0 1]); yticks([0,1]);

load T_cones_sp T_cones_sp S_cones_sp
load T_rods T_rods S_rods
load T_melanopsin T_melanopsin S_melanopsin

plot(SToWls(S_cones_sp),T_cones_sp)
plot(SToWls(S_rods),T_rods)
plot(SToWls(S_melanopsin),T_melanopsin)

ylabel('Sensitivity')

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

% Find peak sensitivities (sure there's a way to do this in less lines)
[~,ms.m_idx(1:3)] = max(T_cones_sp'); %ms = max sensitivities
[~,ms.m_idx(4)]   = max(T_rods'); 
[~,ms.m_idx(5)]   = max(T_melanopsin'); 

ms.fS_sp  = SToWls(S_cones_sp); %full S
ms.fS_r   = SToWls(S_rods);
ms.fS_mel = SToWls(S_melanopsin);

ms.m(1:3) = ms.fS_sp(ms.m_idx(1:3)); %max
ms.m(4)   = ms.fS_r(ms.m_idx(4));
ms.m(5)   = ms.fS_mel(ms.m_idx(5));

%% Show the effect of changing the weight of PC2
figure(1)
cols = lines(3);
xlim([380 780]); xticks(380:100:780);
ylim([-0.3 0.3]); yticks([-0.3,0,0.3]);
ylabel('Power');

cla
fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
    'k','LineStyle','none','FaceAlpha','0.03');

plot(SToWls(S_SPD),pc_p.coeff(:,2),'Color',cols(2,:))
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

plot([ms.m(3),ms.m(3)],[min(ylim),max(ylim)],'Color','k')
plot([ms.m(5),ms.m(5)],[min(ylim),max(ylim)],'Color','k')
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

cla
fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
    'k','LineStyle','none','FaceAlpha','0.03');
plot(SToWls(S_SPD),pc_p.coeff(:,2)*-1,'Color',cols(2,:)*((-1+2)/4))
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

plot([ms.m(3),ms.m(3)],[min(ylim),max(ylim)],'Color','k')
plot([ms.m(5),ms.m(5)],[min(ylim),max(ylim)],'Color','k')
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

cla
fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
    'k','LineStyle','none','FaceAlpha','0.03');
for i = -2:0.5:2
    plot(SToWls(S_SPD),pc_p.coeff(:,2)*i,'Color',cols(2,:)*((i+2)/4))
end
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%%

load sur_vrhel.mat

figure('Position',[plot_where 800 800]), hold on
set(gca, 'FontSize', 16)
set(gcf,'defaultLineLineWidth',2)

%sur_vrhel_n = sur_vrhel(:,[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141]);
sur_vrhel_n = sur_vrhel(:,[1:44,65,69,81:154]);
sur_vrhel_n_c = corr(sur_vrhel_n');
surf(sur_vrhel_n_c*100,'EdgeColor','none')
axis image
colormap gray
S_refs_f = SToWls(S_vrhel);
xticks(6:50:200)
yticks(6:50:200)
set(gca,'XTickLabel',S_refs_f(xticks))
set(gca,'YTickLabel',S_refs_f(yticks))
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

plot3([26,26],[min(ylim),max(ylim)],[max(zlim),max(zlim)],'Color','b')
plot3([50,50],[min(ylim),max(ylim)],[max(zlim),max(zlim)],'Color','k')
plot3([min(xlim),max(xlim)],[26,26],[max(zlim),max(zlim)],'Color','b')
plot3([min(xlim),max(xlim)],[50,50],[max(zlim),max(zlim)],'Color','k')
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

plot3([88,88],[min(ylim),max(ylim)],[max(zlim),max(zlim)],'Color','r')
plot3([78,78],[min(ylim),max(ylim)],[max(zlim),max(zlim)],'Color','g')
plot3([58,58],[min(ylim),max(ylim)],[max(zlim),max(zlim)],'Color',[0.5,0.5,0.5])
plot3([min(xlim),max(xlim)],[88,88],[max(zlim),max(zlim)],'Color','r')
plot3([min(xlim),max(xlim)],[78,78],[max(zlim),max(zlim)],'Color','g')
plot3([min(xlim),max(xlim)],[58,58],[max(zlim),max(zlim)],'Color',[0.5,0.5,0.5])
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% Plot normalised spectral reflectances
figure(1)

cla
hold on
xlim([390 730]); xticks(400:100:700);
ylim([0 10]); yticks([0,10]);
xlabel('Wavelength'); ylabel('Normalised  Reflectance');

plot([440,440],[min(ylim),max(ylim)],'b:')
plot([488,488],[min(ylim),max(ylim)],'k:')

plot(SToWls(S_SRF),T_SRF./T_SRF(26,:))
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure


%%

% figure(2)
% cla
% 
% scatter(T_SRF(26,:),T_SRF(76,:))
% xlim('auto')
% ylim('auto')
% xticks('auto')
% yticks('auto')

%%
[cs,pc_p,plt_lbls] = melcomp_3;

% L against PC1
figure(1)
cla
scatter(squeeze(mean(cs(1,:,:),2)),pc_p.SCORE(:,1),'k.')
xlabel(plt_lbls{1})
ylabel('PC 1')
xlim('auto')
xticks('auto')
ylim('auto')
yticks('auto')

f_SI = polyfit(squeeze(mean(cs(1,:,:),2)),pc_p.SCORE(:,1),1);
x = linspace(min(squeeze(mean(cs(1,:,:),2))),max(squeeze(mean(cs(1,:,:),2))));
y = x*f_SI(1)+f_SI(2);
plot(x,y)

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

% S/I against PC2
cla
scatter(squeeze(mean(cs(18,:,:),2)),pc_p.SCORE(:,2),'k.')
xlabel(plt_lbls{18})
ylabel('PC 2')
xlim('auto')
xticks('auto')
ylim('auto')
yticks('auto')

f_SI = polyfit(squeeze(mean(cs(18,:,:),2)),pc_p.SCORE(:,2),1);
x = linspace(min(squeeze(mean(cs(18,:,:),2))),max(squeeze(mean(cs(18,:,:),2))));
y = x*f_SI(1)+f_SI(2);
plot(x,y)

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% 3D plot - basic (same as previous but without line and with a third dimension)
cla
scatter3(squeeze(mean(cs(18,:,:),2)),pc_p.SCORE(:,2),pc_p.SCORE(:,1),'k.')
%set(gca,'Color',repmat(c_norm(18,2),3,1))
xlabel(plt_lbls{18})
ylabel('PC 2')
zlabel('PC 1')

%% 3D plot - log(S/I) against sqrt(normalised(PC1))

SoI = squeeze(mean(cs(18,:,:),2)); %signal of interest, OR S-over-I (take your pick!)
SoIL = log2(SoI); %Log SoI
pc1 = pc_p.SCORE(:,1);
pc2 = pc_p.SCORE(:,2);
pc3 = pc_p.SCORE(:,3);
pc1n = (pc1-min(pc1)).^(1/4); %pc1 normalised

cla
scatter3(SoIL,pc2,pc1n,'b.')
xlabel('log2(S/I)')
ylabel('PC2')
zlabel('(PC1-min(PC1)).^(1/4)','Interpreter','none')

%% Plot showing segmentation by PC1

plt_seg = 1;

NoD = 30; %nnumber of divisions
m1 = max(pc1n);
cols = lines(NoD); %cols = jet(n);

sc_t = [SoIL pc2]; %scatter temp
sc = NaN([size(sc_t) NoD]);

if plt_seg
    figure('Position',[plot_where 800 800]), hold on
    set(gca, 'FontSize', 16)
    set(gcf,'defaultLineLineWidth',2)
    
    zlabel('(PC1-min(PC1)).^(1/4)','Interpreter','none')
    axis equal
    xlabel('log2(S/I)')
    ylabel('PC2')
    cla
    
    for i=1:NoD
        block(i,[1 2]) = [(i-1)*(m1/NoD), i*(m1/NoD)];
        sc(and(pc1n>=(block(i,1)),pc1n<block(i,2)),:,i) = sc_t(and(pc1n>=(block(i,1)),pc1n<block(i,2)),:);
        fit_t(i,:) = polyfit(sc(~isnan(sc(:,1,i)),1,i),sc(~isnan(sc(:,2,i)),2,i),1); %fitty, lol. Seriously though, temp.
        
        x_seg = linspace(min(sc(:,1,i)),max(sc(:,1,i)),20);
        y = (fit_t(i,1) * x_seg) + fit_t(i,2);
        
        if any(fit_t(i,:))
            scatter3(sc(:,1,i),sc(:,2,i),pc1n,[],cols(i,:));
            plot3(x_seg,y,repmat(mean(block(i,:)),size(x_seg,2),1),'Color',cols(i,:));
        end
        drawnow
        %pause(0.2)
    end
end

%%

plt_ass = 1;

if plt_ass
    figure('Position',[plot_where 800 800]), hold on
    set(gca, 'FontSize', 16)
    set(gcf,'defaultLineLineWidth',2)
    
    x_ass = mean(block(and(fit_t(:,1),fit_t(:,2)),:),2);
    y_1 = fit_t(and(fit_t(:,1),fit_t(:,2)),1);
    y_2 = fit_t(and(fit_t(:,1),fit_t(:,2)),2);
    plot(x_ass,y_1,'r')
    plot(x_ass,y_2,'b')
    
    p_1 = polyfit(log(x_ass),log(y_1),2);
    p_2 = real(polyfit(log(x_ass),log(y_2),2)); %goes imaginary due to negative numbers (I think?)
    y_1p = polyval(p_1, log(x_ass));
    y_2p = polyval(p_2, log(x_ass));
    
    plot(x_ass,exp(y_1p),'r:')
    plot(x_ass,exp(y_2p),'b:')
    
    
    legend('m','c',num2str(p_1),num2str(p_2))
end

%%

return

%%

figure(1)

scatter3(log2(squeeze(mean(cs(13,:,:),2))),pc2,pc1n,'r.')
scatter3(log2(squeeze(mean(cs(14,:,:),2))),pc2,pc1n,'g.')

xlabel(' ')

legend('log(S/I)','log(L/M)','log(L/S)')

%%

close all

for i=[6:9,13:21]
    figure,
    scatter3(log2(squeeze(mean(cs(i,:,:),2))),pc2,pc1n)
    title(plt_lbls{i})
    axis equal
    
    xlabel('target signal')
    ylabel('PC2')
    zlabel('(PC1-min(PC1)).^(1/4)','Interpreter','none')
    view(2)
end

%%

close all

for i=[6:9,13:21]
    figure,
    scatter3(log2(squeeze(mean(cs(i,:,:),2))),pc3,pc1n)
    title(plt_lbls{i})
    axis equal
    
    xlabel('target signal')
    ylabel('PC3')
    zlabel('(PC1-min(PC1)).^(1/4)','Interpreter','none')
    view(2)
end

%%