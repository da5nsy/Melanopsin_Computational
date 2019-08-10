clear, clc, close all

%% Pre-flight

% Display Settings
plt.disp = 1;         % Display figures?
d.s=25;               % display, size
d.MFA = 0.2;          % Marker Face Alpha
d.mktrns = 0.3;       % Marker transparency
set(groot,'defaultfigureposition',[100 100 500 400])
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesFontName', 'Courier')
set(groot,'defaultAxesFontSize',12)
set(groot,'defaultFigureRenderer', 'painters') %renders pdfs as vector graphics
set(groot,'defaultfigurecolor','white')
cols = hsv(10); rng(2);
set(groot,'defaultAxesColorOrder',cols(randperm(size(cols,1)),:))

% MB settings
min_l_scale = 0.62;
max_l_scale = 0.82;
max_s_scale = 0.05;

% Plot saving settings
plt.print = 1; % Save figures
base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\predictingChromaticity';

%% Load Data

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_1',...
    'SSF','SS10',...
    'lum','CIE_10',...
    'mel_offset',0);

%% Colorimetry

[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

% generate colorimetry for the light source itself
[LMSRI_neutral, lsri_neutral] = melcomp_colorimetry(T_SPD, ones(S_sh(3),1), T_SSF, T_lum, S_sh);

%% Basic MB plot

if plt.disp
    figure, hold on
    drawChromaticity('MB10')
    
    % Plot surfaces
    dMECalt = hsv(size(T_SRF,2));    
    for i=1:size(T_SRF,2)
        scatter(lsri(1,i,:),lsri(2,i,:),d.s,'filled','MarkerFaceAlpha',d.MFA)
    end
    
    % Plot illuminant only
    scatter(lsri_neutral(1,:),lsri_neutral(2,:),d.s,'k','filled','MarkerFaceAlpha',d.MFA)
    

    xticks([min(xlim),max(xlim)])
    yticks([min(ylim),max(ylim)])
end

if plt.print
    save2pdf([base,'\BasicMB_1.pdf'])
end

xlim([min_l_scale max_l_scale]);
ylim([0 max_s_scale]);
cleanTicks

if plt.print
    save2pdf([base,'\BasicMB_2.pdf'])
end

%%










