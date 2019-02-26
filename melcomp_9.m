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
    'SRF','Vrhel_nat_2',...
    'SSF','SS10',...
    'mel_offset',0,...
    'lum','CIE_10');
refs=[38, 15, 134, 137, 138, 65, 19, 24, 140, 26];

%% Compute colorimetry of reflectance samples

[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

%compute colours for display
pltc_alt = repmat(jet(size(T_SRF,2))',1,1,size(T_SPD,2));
rng(7); pltc_alt=pltc_alt(:,randperm(size(T_SRF,2)),:); %best combo for differentiating close chromaticities

figure('Position',[plot_where plot_size]), hold on
scatter(lsri(1,:),lsri(2,:),...
    [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)

%rescale diagram
xlim([min_l_scale max_l_scale])
ylim([0 max_s_scale])
xticks([min_l_scale max_l_scale])
yticks([0 max_s_scale])

%% Normalise data for std and set mean to 0

lsri_c = lsri(:,:); %would be nice to not do this transformation, for easier access later, but I'd need to be very careful to check that it didn't change the function of the calculation below

for i=1:4
    lsri_c(i,:) = (lsri_c(i,:) - mean(lsri_c(i,:)))./std(lsri_c(i,:));
end

figure, 
scatter(lsri_c(1,:),lsri_c(2,:),...
    [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
axis equal

lsri_c = reshape(lsri_c,size(lsri));

%% Calculate corrections

% Grey world

lsri_m2 = mean(lsri_c,2);
% hold on
% scatter(lsri_m2(1,:),lsri_m2(2,:),'k.')

lsri_gw = lsri_c - repmat(lsri_m2,1,10,1);
% figure, hold on
% scatter(lsri_gw(1,:),lsri_gw(2,:),...
%     [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
% scatter(0,0,'k.')

% Bright is white

whites =  max(lsri_c,[],2);
lsri_bw = lsri_c - repmat(whites,1,10,1);
figure, hold on
scatter(lsri_bw(1,:),lsri_bw(2,:),...
    [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
scatter(0,0,'k.')


% Melanopsin

%% Compute 'success'



