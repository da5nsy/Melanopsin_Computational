function [MB1_minSD,MB2_minSD,fac1,fac2] = transformToIllIndSpace(offset,wholeset,disp,print)


%% Pre-flight

try
    nargin;
catch
    clear, clc, close all;
    offset = 0;
    wholeset = 0;
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
    plt.print = 1; % Save
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

%% Finding k to minimise SD of entire set

MBx=lsri;
MBx1std=[]; %Initialise variables so that they can be added to
MBx2std=[];

if wholeset %minimise SD for the whole set
    for S1= -2:0.01:2
        MBx(1,:,:)=lsri(1,:,:)+S1*lsri(4,:,:);
        MBx1std=[MBx1std,[S1;std(MBx(1,:))]];
    end
    
    for S2= -2:0.01:2
        MBx(2,:,:)=lsri(2,:,:)+S2*lsri(4,:,:);
        MBx2std=[MBx2std,[S2;std(MBx(2,:))]];
    end
else %minimise SD for each set
    for S1= -2:0.01:2
        MBx(1,:,:)=lsri(1,:,:)+S1*lsri(4,:,:);
        MBx1std=[MBx1std,[S1;mean(std(squeeze(MBx(1,:,:))'))]];
    end
    
    for S2= -2:0.01:2
        MBx(2,:,:)=lsri(2,:,:)+S2*lsri(4,:,:);
        MBx2std=[MBx2std,[S2;mean(std(squeeze(MBx(2,:,:))'))]];
    end
end

if plt.disp
    figure, hold on
    plot(MBx1std(1,:),MBx1std(2,:))
    plot(MBx2std(1,:),MBx2std(2,:))
    xlabel('k')
    ylabel('SD')
    legend('k1','k2')
end

if plt.print
    if wholeset
        save2pdf([base,'\kvsSD.pdf'])
    else
        save2pdf([base,'\kvsSD2.pdf'])
    end
end

% % Apply factors

[MB1_minSD, MB1_minSD_loc] = min(MBx1std(2,:));
[MB2_minSD, MB2_minSD_loc] = min(MBx2std(2,:));

fac1 = MBx1std(1,MB1_minSD_loc);
fac2 = MBx2std(1,MB2_minSD_loc);

MB_star = zeros(size(lsri(1:2,:,:)));
MB_star(1,:,:)=lsri(1,:,:) + fac1 * lsri(4,:,:);
MB_star(2,:,:)=lsri(2,:,:) + fac2 * lsri(4,:,:);

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
