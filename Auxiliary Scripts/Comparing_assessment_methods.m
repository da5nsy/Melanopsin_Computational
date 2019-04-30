% Different types of transform

%% Pre-flight
clear, clc, close all

plot_where = [500,200];
plot_size  = [800,400];

min_l_scale = 0.62;
max_l_scale = 0.82;
max_s_scale = 0.04;

set(0,'defaultAxesFontName', 'Courier')

%% Load Data

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_1',...
    'SSF','SS10',...
    'mel_offset',0,...
    'lum','CIE_10');
refs=[38, 15, 134, 137, 138, 65, 19, 24, 140, 26];

%% Compute colorimetry of reflectance samples

[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

%compute colours for display
pltc_alt = repmat(jet(size(T_SRF,2))',1,1,size(T_SPD,2));
rng(7); pltc_alt=pltc_alt(:,randperm(size(T_SRF,2)),:); %best combo for differentiating close chromaticities

% figure('Position',[plot_where plot_size]), hold on
% scatter(lsri(1,:),lsri(2,:),...
%     [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)

% %rescale diagram
% xlim([min_l_scale max_l_scale])
% ylim([0 max_s_scale])
% xticks([min_l_scale max_l_scale])
% yticks([0 max_s_scale])

%% Normalise data for std and set mean to 0

% % Testing standardness of distributions
% hk=10;
% figure,hist(lsri_c(1,:),hk)
% figure,hist(lsri_c(2,:),hk)
% 
% figure,hist(log(lsri_c(1,:)),hk)
% figure,hist(log(lsri_c(2,:)),hk)

lsri = log(lsri);

lsri_c = lsri(:,:); %would be nice to not do this transformation, for easier access later, but I'd need to be very careful to check that it didn't change the function of the calculation below

for i=1:size(lsri,1)
    lsri_c(i,:) = (lsri_c(i,:) - mean(lsri_c(i,:)))./std(lsri_c(i,:));
end

% figure, 
% scatter(lsri_c(1,:),lsri_c(2,:),...
%     [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
% axis equal

lsri_c = reshape(lsri_c,size(lsri));




%% Grey world

lsri_m2 = mean(lsri_c,2);
% hold on
% scatter(lsri_m2(1,:),lsri_m2(2,:),'k.')

lsri_gw = lsri_c - repmat(lsri_m2,1,size(lsri,2),1);
figure, hold on
scatter(lsri_gw(1,:),lsri_gw(2,:),...
    [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
scatter(0,0,'k.')

axis equal
xlim([-3 3])
ylim([-3 3])
cleanTicks

title('Grey world correction')

%% Melanopsin

[sf_l,sf_s] = melcomp_6_calcsf(lsri_c,0:0.01:1,-2:0.01:-0.5,1,pltc_alt); %calculates scaling factors

lsri_mel = [lsri_c(1,:,:)+sf_l*lsri_c(4,:,:);lsri_c(2,:,:)+sf_s*lsri_c(4,:,:)];

figure
scatter(lsri_mel(1,:),lsri_mel(2,:),...
    [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)

axis equal
xlim([-3 3])
ylim([-3 3])
cleanTicks

title('Melanopsin-based correction')

%% Forsyth measure

ForsythMeasurementOfSuccess(lsri_gw,pltc_alt)
ForsythMeasurementOfSuccess(lsri_mel,pltc_alt)

%% kmeans

KMeansMark(lsri_gw)
KMeansMark(lsri_mel)


