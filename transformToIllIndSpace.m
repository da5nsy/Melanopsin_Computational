function [minSD_l,minSD_s,sf_l,sf_s] = transformToIllIndSpace(offset,wholeset,disp,print)


%% Pre-flight

try
    nargin;
catch
    clear, clc, close all;
    offset = 75;
    wholeset = 1;
end

% Display Settings
if exist('disp','var') % Display figures?
    plt.disp = disp;
else
    plt.disp = 1; 
end         
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
if exist('print','var')
    plt.print = print;
else
    plt.print = 0; % Save
end
if plt.print
    warning('plt.print is enabled - you sure? This will overwrite existing figures.')
end
base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\transformToIllIndSpace';

%% Load Data

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_1',...
    'SSF','SS10',...
    'lum','CIE_10',...
    'mel_offset',offset);

%% Colorimetry

[~, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

% Log tranform, zero mean, and normalise by SD
lsri = log(lsri);
lsri = lsri - mean(lsri(:,:),2);
lsri = lsri./std(lsri(:,:),[],2);

%% Finding k to minimise SD of entire set

if plt.disp
    plt.sfs = 1;
end

if wholeset
    [sf_l,sf_s,minSD_l,minSD_s] = calcsf(lsri, -2:0.01:2, -2:0.01:2,plt,1);
else
    [sf_l,sf_s,minSD_l,minSD_s] = calcsf(lsri, -2:0.01:2, -2:0.01:2,plt,0);
end

if plt.print
    if wholeset
        save2pdf([base,'\kvsSD.pdf'])
    else
        save2pdf([base,'\kvsSD2.pdf'])
    end
end

% % Apply factors

MB_star = zeros(size(lsri(1:2,:,:)));
MB_star(1,:,:)=lsri(1,:,:) + sf_l * lsri(4,:,:);
MB_star(2,:,:)=lsri(2,:,:) + sf_s * lsri(4,:,:);

if plt.disp    
    figure, hold on    
    for i=1:size(T_SRF,2)
        scatter(MB_star(1,i,:),MB_star(2,i,:),d.s,'filled','MarkerFaceAlpha',d.MFA)
    end
    cleanTicks
    xlabel('{\itl}_{MB} + {\itk_1i}_{MB}');
    ylabel('{\its}_{MB} + {\itk_2i}_{MB}');
end

if plt.print
    if wholeset
        save2pdf([base,'\correctedChromaticities.pdf'])
    else
        save2pdf([base,'\correctedChromaticities2.pdf'])
    end
end

end
