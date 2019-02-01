%% Creating visual narrative for melcomp presentation

% I gave a presentation on 2019/01/29 to Hannah Smithson's group in Oxford
% on the latest progress of my computational work. Here is where I
% generated the figures for that presentation.

%Pre-flight
clear, clc, close all

plot_where = [500,200];
plot_size  = [800,400];

min_l_scale = 0.62;
max_l_scale = 0.82;
max_s_scale = 0.04;

set(0,'defaultAxesFontName', 'Courier')

base = 'C:\Users\cege-user\Dropbox\UCL\Presentations\20190129 Computational Study - Oxford\figs';
print_figures = 1;

%% Load Data

% Load observer data

load T_cones_ss10.mat;
T_SSF = T_cones_ss10;
S_SSF = S_cones_ss10;
clear T_cones_ss10 S_cones_ss10

load T_rods T_rods S_rods
load T_melanopsin T_melanopsin S_melanopsin
T_mel = SplineCmf(S_melanopsin, T_melanopsin, S_melanopsin - [10, 0, 0],1); %Increasing the range of this function in case it ends up limiting the range of S_sh, and shorten variable names
S_mel = S_melanopsin - [10, 0, 0];

% Load reflectance data

load sur_vrhel
refs=[38, 15, 134, 137, 138, 65, 19, 24, 140, 26];
sur_vrhel_n = sur_vrhel(:,refs);
T_SRF = sur_vrhel_n;
S_SRF = S_vrhel;

% Load daylight data

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)
T_SPD=final; clear final
T_SPD = T_SPD(:,1:20:end);
S_SPD=[300,5,161];

%% reduce all data down to the common range/interval

S_sh = [max([S_SPD(1),S_SRF(1),S_SSF(1),S_rods(1),S_mel(1)]),...
    max([S_SPD(2),S_SRF(2),S_SSF(2),S_rods(2),S_mel(2)]),...
    min([S_SPD(3),S_SRF(3),S_SSF(3),S_rods(3),S_mel(3)])];
%S_shared: work out what the lowest common denominator for the range/interval of the data is

T_SPD = SplineSpd(S_SPD,T_SPD,S_sh);
T_SRF = SplineSrf(S_SRF,T_SRF,S_sh,1); %ended with same value
T_SSF  = SplineCmf(S_SSF,T_SSF,S_sh)';
T_rods = SplineCmf(S_rods,T_rods,S_sh)';
T_mel  = SplineCmf(S_mel,T_mel,S_sh)';
[S_SPD, S_SRF, S_SSF, S_rods, S_mel] = deal(S_sh);

% combine sensitivity vectors
T_LMSRI=[T_SSF,T_rods,T_mel];
S_LMSRI=S_sh;

%% Plot MB chromaticity diagram

% compute chromaticities of points on spectral locus
sf_10 = [0.69283932, 0.34967567, 0.05547858]; %energy 10deg from CIE 170-2:2015
T_lum = sf_10(1)*T_SSF(:,1)+sf_10(2)*T_SSF(:,2);
spectral_locus = LMSToMacBoynDG(T_SSF',T_SSF',T_lum');

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
xlabel('{\itl}_{MB}');
ylabel('{\its}_{MB}');

if print_figures
    save2pdf([base,'\MB.pdf'])
end

%black version
scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
if print_figures
    save2pdf([base,'\MBblack.pdf'])
end

%% Compute colorimetry of reflectance samples

T_rad = zeros([S_sh(3),size(T_SRF,2),size(T_SPD,2)]);
LMSRI = zeros([size(T_LMSRI,2),size(T_SRF,2),size(T_SPD,2)]);
lsri  = zeros([4,size(T_SRF,2),size(T_SPD,2)]);
t_r   = zeros([2,size(T_SRF,2),size(T_SPD,2)]); %t for temp
t_i   = zeros([2,size(T_SRF,2),size(T_SPD,2)]); %t for temp

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_SRF.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_LMSRI'*T_rad(:,:,i);
    lsri(1:2,:,i) = LMSToMacBoynDG(LMSRI(1:3,:,i),T_SSF',T_lum');
    t_r(:,:,i)    = LMSToMacBoynDG(LMSRI([1,2,4],:,i),[T_SSF(:,1:2)';T_rods'],T_lum');
    t_i(:,:,i)    = LMSToMacBoynDG(LMSRI([1,2,5],:,i),[T_SSF(:,1:2)';T_mel'],T_lum');
end
lsri(3,:,:) = t_r(2,:,:); clear t_r
lsri(4,:,:) = t_i(2,:,:); clear t_i

%compute colours for display
pltc_alt = repmat(jet(size(T_SRF,2))',1,1,size(T_SPD,2));
rng(7); pltc_alt=pltc_alt(:,randperm(size(T_SRF,2)),:); %best combo for differentiating close chromaticities

%plot MB with points, not normalised
rng(1); n_ill = randi(size(T_SPD,2)); %pick a random spectrum (55th under current settings, changes if you downsample the illuminant data, for example)
scatter(lsri(1,:,n_ill),lsri(2,:,n_ill),[],pltc_alt(:,:,n_ill)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
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
labels_vrhel(137).label = 'peach skin -- yellow';
for i=1:length(refs)
    text(lsri(1,i,n_ill)+0.005,lsri(2,i,n_ill)+0.00015,labels_vrhel(refs(i)).label,'Rotation',10,'FontName','Courier')
end
if print_figures
    save2pdf([base,'\MBsingleset_n_text.pdf'])
end

%scatter(0.699237,0.025841,'rs') EE white

%% Plot chromaticities of reflectance samples under all illums

% all points, in colour
cla
scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
scatter(lsri(1,:),lsri(2,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
if print_figures
    save2pdf([base,'\MBall.pdf'])
end

%all points in grey
cla
scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
scatter(lsri(1,:),lsri(2,:),[],[0.5,0.5,0.5],'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
if print_figures
    save2pdf([base,'\MBallgrey.pdf'])
end

%with means under each illuminant
lsri_m2 = mean(lsri,2);
scatter(squeeze(lsri_m2(1,:,:)),squeeze(lsri_m2(2,:,:)),'r*')
if print_figures
    save2pdf([base,'\MBmeans.pdf'])
end

%% Grey world adaptation
corrector = lsri_m2-lsri_m2(:,:,1);
lsri_mc = lsri_m2 - corrector;
lsri_c = lsri - corrector;
cla
scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
scatter(lsri_c(1,:),lsri_c(2,:),[],[0.5,0.5,0.5],'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
scatter(lsri_mc(1,1,1),lsri_mc(2,1,1),'r*')
xlabel('{\it }_{ }');
ylabel('{\it }_{ }');
if print_figures
    save2pdf([base,'\MBc.pdf'])
end

%% Splits

figure('Name','split','Position',[plot_where plot_size]), hold on

s(1) = subplot(1,2,1);
hold on
s(2) = subplot(1,2,2);
hold on

xlim(s(1),[min_l_scale max_l_scale])
xticks(s(1),[min_l_scale max_l_scale])
xlabel(s(1),'{\itl}_{MB}');
yticks(s(1),[])

xlim(s(2),[0 max_s_scale])
xticks(s(2),[0 max_s_scale])
xlabel(s(2),'{\its}_{MB}');
yticks(s(2),[])

if print_figures
    save2pdf([base,'\split_empty.pdf']) %blank one for clarity of introduction in presentation
end

ylabel(s(1),'{\itI}');
scatter(s(1),lsri(1,:),LMSRI(5,:),[],[0.5,0.5,0.5],'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
scatter(s(2),lsri(2,:),LMSRI(5,:),[],[0.5,0.5,0.5],'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
yticks(s(1),[min(ylim),max(ylim)])
yticks(s(2),[])

save2pdf([base,'\split_I_grey.pdf'])

%LMSRI (with category colours)
plotOrderNums = [5,1,2,3];
plotOrderNames = {'I','L','M','S'};
for i=1:length(plotOrderNums)
    cla(s(1))
    cla(s(2))
    
    ylabel(s(1),['{\it',plotOrderNames{i},'}']);
    
    scatter(s(1),lsri(1,:),LMSRI(plotOrderNums(i),:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    scatter(s(2),lsri(2,:),LMSRI(plotOrderNums(i),:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    
    yticks(s(1),[min(ylim),max(ylim)])
    yticks(s(2),[])
    if print_figures
        save2pdf([base,'/split_',plotOrderNames{i},'.pdf'])
    end
end

%lsri
plotOrderNums = [1,2,4];
plotOrderNames = {'{\itl}_{MB}','{\its}_{MB}','{\iti}_{MB}'};
plotOrderSaveNames = {'lmb','smb','imb'};
for i=1:length(plotOrderNums)
    cla(s(1))
    cla(s(2))
    
    ylabel(s(1),plotOrderNames{i});
    
    scatter(s(1),lsri(1,:),lsri(plotOrderNums(i),:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    scatter(s(2),lsri(2,:),lsri(plotOrderNums(i),:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    
    yticks(s(1),[min(ylim),max(ylim)])
    yticks(s(2),[])
    
    if print_figures
        save2pdf([base,'/split_',plotOrderSaveNames{i},'.pdf'])
    end
end

%% Calibartion by l, s or i

figure(1)

cla
scatter(lsri(1,:)-lsri(1,:),lsri(2,:)-lsri(1,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
axis('auto')
xticks([])
yticks([])

xlabel('{\itl}_{MB} - {\itl}_{MB}');
ylabel('{\its}_{MB} - {\itl}_{MB}');

if print_figures
    save2pdf([base,'\MBminMB1.pdf'])
end

cla
scatter(lsri(1,:)-lsri(2,:),lsri(2,:)-lsri(2,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
axis('auto')
xticks([])
yticks([])

xlabel('{\itl}_{MB} - {\its}_{MB}');
ylabel('{\its}_{MB} - {\its}_{MB}');

if print_figures
    save2pdf([base,'\MBminMB2.pdf'])
end

[sf_l,sf_s] = melcomp_6_calcsf(lsri); %calculates scaling factors
cla
scatter(lsri(1,:)+sf_l*lsri(4,:),lsri(2,:)+sf_s*lsri(4,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
axis('auto')
xticks([])
yticks([])

xlabel('{\itl}_{MB} - {\itk_1i}_{MB}');
ylabel('{\its}_{MB} - {\itk_2i}_{MB}');

if print_figures
    save2pdf([base,'\MBminMB3.pdf'])
end

%% Plot spectral reflectance functions

figure('Position',[plot_where plot_size],'defaultLineLineWidth',4), hold on
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
plot(SToWls(S_SPD),pc_p.coeff(:,1:3))
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
    pc(i) = melcomp_6_looper(pca_range(i),1,0);
    disp(pca_range(i))
end

[~, mel_peak_loc] = max(T_melanopsin);
S_melanopsin_f = SToWls(S_melanopsin);
mel_peak = S_melanopsin_f(mel_peak_loc);

for i=1:length(pca_range)
    ex(i) = pc(i).explained(3);
end

%figure, hold on
plot(pca_range+mel_peak,ex/max(ex),'k','DisplayName','PC3 score')
%plot

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
plot(SToWls(S_SSF),T_SSF)
plot(SToWls(S_melanopsin),T_melanopsin)

if print_figures
    save2pdf([base,'\optwithSFF.pdf'])
end

%% Plot 3D at first peak

[~,pks_locs] = findpeaks(ex);

% for i = 1:length(pks_locs)
%     melcomp_6_looper(pca_range(pks_locs(i)),1,1)
% end

figure('Position',[plot_where plot_size],'defaultLineLineWidth',4), hold on
xlabel('{\itl}_{MB}');
ylabel('{\its}_{MB}');
zlabel('{\iti}_{MB}');
melcomp_6_looper(pca_range(pks_locs(1)),1,1)

view(3)

xticks([min(xlim),max(xlim)]);
yticks([min(ylim),max(ylim)]);
zticks([min(zlim),max(zlim)]);

if print_figures
    save2pdf([base,'\3D.pdf'])
end

%% Show impact of varying parameters

d = dir('Auxiliary Scripts\Optimality variation\opt*.mat');

figure(3), hold on
cla
for i = 1:length(d)
    clear out
    load([d(i).folder,'\', d(i).name]);
    plot(out(1,:),out(2,:),'-','DisplayName',d(i).name);
end

legend('off')
ylim('auto')

if print_figures
    save2pdf([base,'\optopt.pdf'])
end

