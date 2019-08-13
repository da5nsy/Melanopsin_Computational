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
plt.print = 0; % Save 
if plt.print
    warning('plt.print is enabled - you sure? This will overwrite existing figures.')
end
base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\comparisonFourAlgos';

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

%% Perform CC algos

Lum = zeros(size(T_SRF,2),size(T_SPD,2));
for i=1:size(lsri,3)
    Lum(:,i) = T_lum'*(T_SRF.*T_SPD(:,i));
end

[output,sf_l,sf_s] = performCC(lsri,Lum,1);

%% Plot outputs

figure('Position',[100 100 800 1000],'Renderer','opengl') %~500kb vs ~150kb with opengl (because it makes it a bitmap rather than a vector

for i = 1:size(output,4)
    subplot(2,2,i)
    hold on
    for j=1:size(T_SRF,2)
        scatter(output(1,j,:,i),output(2,j,:,i),d.s,'filled','MarkerFaceAlpha',d.MFA)
    end
    cleanTicks
    xlabel('{\itl}_{MB}*');
    ylabel('{\its}_{MB}*');
end

if plt.print
    save2pdf([base,'\output2.pdf'])
end

%% Kmeans

for i = 1:size(output,4)
    KMeansMark(squeeze(output(:,:,:,i)))
end

%% Perform corrections

% Calculate luminance for BiW
Lum = zeros(size(T_SRF,2),size(T_SPD,2));
for i=1:size(lsri,3)
    Lum(:,i) = T_lum'*(T_SRF.*T_SPD(:,i));
end

pcSurfRange = 20:10:100;

rng(1)


%%
figure, plot(pcSurfRange,mark','o')
legend('Location','best')
ylim([0 1])







