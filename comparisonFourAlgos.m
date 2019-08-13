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
    'SRF','Vrhel_nat_1',... %Vrhel_nat_1 or Vrhel_nat_extended_exclSimilar
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

for pcSurf = [100,50]
    [output,mark,sf_l,sf_s,sel_store] = performCCandMarkForRange(lsri,Lum,0,pcSurf);
    
    % Plot outputs
    
    figure('Position',[100 100 800 1000],'Renderer','opengl') %~500kb vs ~150kb with opengl (because it makes it a bitmap rather than a vector
    
    sel_store = squeeze(sel_store);
    
    for i = 1:size(output,4)
        subplot(2,2,i)
        hold on
        for j=1:size(T_SRF,2)
            pointsToPlot_l = squeeze(output(1,:,:,i));
            pointsToPlot_s = squeeze(output(2,:,:,i));
            scatter(pointsToPlot_l(sel_store == j),pointsToPlot_s(sel_store == j),...
                d.s,'filled','MarkerFaceAlpha',d.MFA)
        end
        cleanTicks
        xlabel('{\itl}_{MB}*');
        ylabel('{\its}_{MB}*');
    end
    
    disp(mark)
    
    if plt.print
        save2pdf([base,'\output',num2str(pcSurf),'.pdf'])
    end
    
end

%% Try a full range

pcSurfRange = 20:10:100;

[output,mark,sf_l,sf_s,sel_store] = performCCandMarkForRange(lsri,Lum,0,pcSurfRange);

figure, plot(pcSurfRange,mark','o')
legend({'DN','GW','BiW','Mel'},'Location','best')
ylim([0 1])

xlabel('% of surfaces')
ylabel('k-means-mark')

if plt.print
    save2pdf([base,'\outputRange.pdf'])
end





