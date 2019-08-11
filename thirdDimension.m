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

% Plot saving settings
plt.print = 1; % Save 
if plt.print
    warning('plt.print is enabled - you sure? This will overwrite existing figures.')
end
base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\thirdDimension';

%% Load Data

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_1',...
    'SSF','SS10',...
    'lum','CIE_10',...
    'mel_offset',0);

%% Colorimetry

[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

% Log tranform, zero mean, and normalise by SD
lsri = log(lsri);
lsri = lsri - mean(lsri(:,:),2);
lsri = lsri./std(lsri(:,:),[],2);

%%

figure, hold on

% Plot surfaces
for i=1:size(T_SRF,2)
    scatter3(lsri(1,i,:),lsri(2,i,:),lsri(4,i,:),d.s,'filled','MarkerFaceAlpha',d.MFA)
end

view(3)

xticks([-2 0 2])
yticks([-2 0 2])
zticks([-2 0 2])

xlabel('{\itl}_{MB}');
ylabel('{\its}_{MB}');
zlabel('{\iti}_{MB}');

if plt.print
    save2pdf([base,'\ZL.pdf'])
end

view(23,-38)

if plt.print
    save2pdf([base,'\viewpoint.pdf'])
end
