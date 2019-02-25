% Playing around with a rubric that marks points based on Forsyth's
% averaging.


%Pre-flight
clear, clc, close all

plot_where = [500,200];
plot_size  = [800,400];

min_l_scale = 0.62;
max_l_scale = 0.82;
max_s_scale = 0.04;

set(0,'defaultAxesFontName', 'Courier')

base = 'C:\Users\cege-user\Dropbox\UCL\Presentations\20190129 Computational Study - Oxford\figs';

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

scatter(lsri(1,:),lsri(2,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)

%rescale diagram
xlim([min_l_scale max_l_scale])
ylim([0 max_s_scale])
xticks([min_l_scale max_l_scale])
yticks([0 max_s_scale])

%% Calculate corrections

% Grey world

% Bright is white

% Melanopsin

%% Compute 'success'



