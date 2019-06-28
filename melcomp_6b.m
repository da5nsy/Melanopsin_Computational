%% Creating visual narrative for melcomp presentation

% I gave a presentation on 2019/01/29 to Hannah Smithson's group in Oxford
% on the latest progress of my computational work. 
%
% This is what I originally wrote the script for.
%
% Following this, I was invited to give a remote seminar to Bevil Conway &
% co at NIH. I have forked the code for this.


%Pre-flight
clear, clc, close all

plot_where = [500,200];
plot_size  = [800,400];

min_l_scale = 0.62;
max_l_scale = 0.82;
max_s_scale = 0.04;

set(0,'defaultAxesFontName', 'Courier')

base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\melcomp_6b_figs';
print_figures = 0;

%% Load Data

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_2',...
    'SSF','SS10',...
    'mel_offset',0,...
    'lum','CIE_10');
refs=[38, 15, 134, 137, 138, 65, 19, 24, 140, 26];

%% Plot MB chromaticity diagram

% compute chromaticities of points on spectral locus
spectral_locus = LMSToMacBoyn(T_SSF(:,1:3)',T_SSF(:,1:3)',T_lum');

% compute display colours for points on spectral locus
load T_xyz1931.mat
T_xyz1931 = SplineCmf(S_xyz1931,T_xyz1931,S_sh);
RGB = XYZToSRGBPrimary(T_xyz1931);
RGB(RGB<0) = 0;
RGB(RGB>1) = 1;

figure('Name','MB','Position',[plot_where plot_size]), hold on
scatter(spectral_locus(1,:),spectral_locus(2,:),[],RGB','filled')
xlim([0.5 1])
ylim([0 1])
xticks([0.5 1])
yticks([0 1])
xlabel('{\itl}_{MB,10}');
ylabel('{\its}_{MB,10}');

if print_figures
    save2pdf([base,'\MB.pdf'])
end

%black version
scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
if print_figures
    save2pdf([base,'\MBblack.pdf'])
end

%% Compute colorimetry of reflectance samples

[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

%compute colours for display
pltc_alt = repmat(jet(size(T_SRF,2))',1,1,size(T_SPD,2));
rng(7); pltc_alt=pltc_alt(:,randperm(size(T_SRF,2)),:); %best combo for differentiating close chromaticities

%plot MB with points, not normalised
rng(1); n_ill = randi(size(T_SPD,2)); %pick a random spectrum (55th under current settings, changes if you downsample the illuminant data, for example)
scatter(lsri(1,:,n_ill),lsri(2,:,n_ill),'k','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
if print_figures
    save2pdf([base,'\MBsingleset_nn.pdf'])
end

%rescale diagram
xlim([min_l_scale max_l_scale])
ylim([0 max_s_scale])
xticks([min_l_scale max_l_scale])
yticks([0 max_s_scale])
if print_figures
    save2pdf([base,'\MBsingleset_n.pdf'])
end

%add labels
load sur_vrhel.mat
labels_vrhel(137).label = 'peach skin -- yellow'; %correct typo
for i=1:length(refs)
    text(lsri(1,i,n_ill)+0.005,lsri(2,i,n_ill)+0.00015,labels_vrhel(refs(i)).label,'Rotation',5,'FontName','Courier')
end
if print_figures
    save2pdf([base,'\MBsingleset_n_text.pdf'])
end

%scatter(0.699237,0.025841,'rs') EE white

%% Plot chromaticities of reflectance samples under all illums

%all points in grey
cla
scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
scatter(lsri(1,:),lsri(2,:),[],[0.5,0.5,0.5],'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
scatter(lsri(1,:,n_ill),lsri(2,:,n_ill),'k','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
for i=1:length(refs)
    text(lsri(1,i,n_ill)+0.005,lsri(2,i,n_ill)+0.00015,labels_vrhel(refs(i)).label,'Rotation',5,'FontName','Courier')
end
if print_figures
    save2pdf([base,'\MBallgrey.pdf'])
end

















%% % - %
%% Plot spectral reflectance functions

figure('Position',[plot_where plot_size],'defaultLineLineWidth',4), hold on
sur_vrhel_n = sur_vrhel(:,refs);
for i=1:length(refs)
    plot(SToWls(S_vrhel),sur_vrhel_n(:,i),'Color',pltc_alt(:,i,1),'DisplayName',char(labels_vrhel(refs(i)).label));
end
xlim([S_vrhel(1), S_vrhel(1)+S_vrhel(2)*S_vrhel(3)])

xticks('auto')
yticks([min(ylim) max(ylim)])

xlabel('Wavelength (nm)')
ylabel('Reflectance (/1)')

l = legend;
l.FontSize = 6;
l.Location = 'northwest';

if print_figures
    save2pdf([base,'\SRF.pdf'])
end

%% Plot SPD pca

[pc_p.coeff, pc_p.score, pc_p.latent, pc_p.tsquared, pc_p.explained, pc_p.mu] = pca(T_SPD');

figure(3), cla
plot(SToWls(S_sh),pc_p.coeff(:,1:3))
xlim([S_sh(1), S_sh(1)+S_sh(2)*S_sh(3)])

xticks('auto')
yticks([min(ylim), 0,  max(ylim)])

xlabel('Wavelength (nm)')
ylabel('Coefficient value')

l = legend({'PC1','PC2','PC3'});
l.FontSize = 6;
l.Location = 'southwest';

if print_figures
    save2pdf([base,'\SPD.pdf'])
end

%% Ref corr square

figure('Position',[plot_where plot_size(1) plot_size(1)])
hold on

vrhel_square = corr(T_SRF');
S_RF_f = SToWls(S_vrhel);

%black out top triangle
for i=1:length(vrhel_square)
    for j=1:length(vrhel_square)
        if i>j
            vrhel_square(i,j) = 0;
        end
    end
end

imagesc(vrhel_square.^7) %power increases visual contrast
axis image
colormap gray
set(gca,'XTickLabel',S_RF_f(xticks))
set(gca,'YTickLabel',S_RF_f(yticks))
colorbar
xlabel('Wavelength (nm)')
ylabel('Wavelength (nm)')

if print_figures
    save2pdf([base,'\VS.pdf'])
end

%% Optimality

figure(3)
cla
clear ex pc %only needed during debugging when rerunning script

pca_range = -70:1:130;

for i=1:length(pca_range)
    pc(i) = melcomp_6_looper('mel_offset',pca_range(i));
    disp(pca_range(i))
end

load T_melanopsin.mat
[~, mel_peak_loc] = max(T_melanopsin);
S_melanopsin_f = SToWls(S_melanopsin);
mel_peak = S_melanopsin_f(mel_peak_loc);

for i=1:length(pca_range)
    ex(i) = pc(i).explained(3);
end

%figure, hold on
plot(pca_range+mel_peak,ex/max(ex),'k','DisplayName','PC3 score')

xlim('auto')
ylim([0 1])
yticks([min(ylim),max(ylim)])
xlabel('Wavelength shift (nm)')
ylabel('PC3 score')
legend('off')

if print_figures
    save2pdf([base,'\opt.pdf'])
end

% add spectral sensitivity functions
plot(SToWls(S_sh),T_SSF);

if print_figures
    save2pdf([base,'\optwithSFF.pdf'])
end

%% Plot 3D at first peak

[~,pks_locs] = findpeaks(ex);

% for i = 1:length(pks_locs)
%     melcomp_6_looper(pca_range(pks_locs(i)),1,1)
% end

figure('Position',[plot_where plot_size],'defaultLineLineWidth',4), hold on
xlabel('{\itl}_{MB,10}');
ylabel('{\its}_{MB,10}');
zlabel('{\iti}_{MB,10}');
melcomp_6_looper('mel_offset',pca_range(pks_locs(1)),'plt',1);

view(3)

xticks([min(xlim),max(xlim)]);
yticks([min(ylim),max(ylim)]);
zticks([min(zlim),max(zlim)]);

if print_figures
    save2pdf([base,'\3D.pdf'])
end

%% Show impact of varying parameters

figure(3), hold on
cla
legend('on')
ylim('auto')
    
pca_range = -70:10:130;

% default
for i=1:length(pca_range)
    pc(i) = melcomp_6_looper('mel_offset',pca_range(i));
    disp(pca_range(i))
end
for i=1:length(pca_range)
    ex(i) = pc(i).explained(3);
end
plot(pca_range,ex,'DisplayName','Standard')

% norm off
for i=1:length(pca_range)
    pc(i) = melcomp_6_looper('mel_offset',pca_range(i),'norm',0);
    disp(pca_range(i))
end
for i=1:length(pca_range)
    ex(i) = pc(i).explained(3);
end
plot(pca_range,ex,'DisplayName','diff MB scaling')

% different refs
for i=1:length(pca_range)
    pc(i) = melcomp_6_looper('mel_offset',pca_range(i),'SRF','Vrhel_nat_1');
    disp(pca_range(i))
end
for i=1:length(pca_range)
    ex(i) = pc(i).explained(3);
end
plot(pca_range,ex,'DisplayName','diff SRF')

% different SSFs
for i=1:length(pca_range)
    pc(i) = melcomp_6_looper('mel_offset',pca_range(i),'SSF','SP','lum','SP');
    disp(pca_range(i))
end
for i=1:length(pca_range)
    ex(i) = pc(i).explained(3);
end
plot(pca_range,ex,'DisplayName','diff SSF')

if print_figures
    save2pdf([base,'\optopt.pdf'])
end


