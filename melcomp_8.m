% Playing around with a rubric that marks points as belonging to a set
% based on euclidian distance.


%Pre-flight
clear, clc, close all

plot_where = [500,200];
plot_size  = [800,400];

min_l_scale = 0.62;
max_l_scale = 0.82;
max_s_scale = 0.04;

set(0,'defaultAxesFontName', 'Courier')

%base = 'C:\Users\cege-user\Dropbox\UCL\Presentations\20190129 Computational Study - Oxford\figs';

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
xlabel('{\itl}_{MB}');
ylabel('{\its}_{MB}');

%% Compute colorimetry of reflectance samples

[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

%compute colours for display
pltc_alt = repmat(jet(size(T_SRF,2))',1,1,size(T_SPD,2));
rng(7); pltc_alt=pltc_alt(:,randperm(size(T_SRF,2)),:); %best combo for differentiating close chromaticities

scatter(lsri(1,:),lsri(2,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)

%rescale diagram
xlim([min_l_scale max_l_scale])
ylim([0 max_s_scale])
xticks([min_l_scale max_l_scale])
yticks([0 max_s_scale])

%% Calculate means

lsri_m = [mean(lsri(1,:,:),3);mean(lsri(2,:,:),3)];
scatter(squeeze(lsri_m(1,:)),squeeze(lsri_m(2,:)),'k*')

%% Calculate distances between each point and each mean.

D = zeros([size(lsri_m,2),size(lsri,2),size(lsri,3)]);

for i=1:size(D,1) % for each mean
    for j=1:size(D,2) %for each sample
        for k=1:size(D,3) %for each illuminant
                D(i,j,k) = sqrt((lsri_m(1,i)-lsri(1,j,k))^2 + (lsri_m(2,i)-lsri(2,j,k))^2);
        end
    end
end

%for each point, list the closest mean

NRST = zeros([size(lsri,2),size(lsri,3)]); %nearest

for i=1:size(lsri,2)
    for j=1:size(lsri,3)
        [~,mloc] = min(D(:,i,j));
        NRST(i,j) = mloc;
    end
end

%figure,hist(NRST'),legend

%%

figure('Name','MB','Position',[plot_where plot_size]), hold on
xlim([min_l_scale max_l_scale])
ylim([0 max_s_scale])
xticks([min_l_scale max_l_scale])
yticks([0 max_s_scale])
xlabel('{\itl}_{MB}');
ylabel('{\its}_{MB}');

scatter(lsri(1,:),lsri(2,:),[],pltc_alt(:,NRST)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)

scatter(squeeze(lsri_m(1,:)),squeeze(lsri_m(2,:)),'k*')

%axis equal

%%

BSLN = repmat([1:10]',1,130); %baseline
CMP = BSLN == NRST; %compare

CLP = mean(CMP(:)); 
%This value indicates the percentage of correctly labelled points
