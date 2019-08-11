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
plt.print = 1; % Save 
if plt.print
    warning('plt.print is enabled - you sure? This will overwrite existing figures.')
end
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

%% First level signals

if plt.disp
    figure, hold on
    drawChromaticity('MB10')
    for i=1:size(T_SRF,2)
        scatter3(lsri_neutral(1,1,:),lsri_neutral(2,1,:),LMSRI(5,i,:),d.s,'filled','MarkerFaceAlpha',d.MFA)
    end
    scatter3(lsri_neutral(1,1,:),lsri_neutral(2,1,:),LMSRI_neutral(5,1,:),d.s,'k','filled','MarkerFaceAlpha',d.MFA)
    
    view(3)
    xlim([min_l_scale, max_l_scale])
    ylim([0, max_s_scale])
    cleanTicks
    xlabel('{\itl}_{MB}');
    ylabel('{\its}_{MB}');
    zlabel('{\itI}_{MB}');
end

if plt.print
    save2pdf([base,'\CvsI.pdf'])
end

%% Second level signals

if plt.disp
    plt.names={'L','M','S','R','I'}; 
    plt.N1 = [1,NaN,NaN,1,2,NaN,1,2,3]; % Plot numbers
    plt.N2 = [2,NaN,NaN,3,3,NaN,5,5,5];
    
    %figure('units','normalized','outerposition',[0 0 1 1])    
    figure('units','normalized','outerposition',[0 0 1 1],'Renderer','opengl')
    
    for i=1:9
        if isnan(plt.N1(i))
            continue
        end
        sp(i)=subplot(3,3,i);
        hold on
        
        for j=1:size(T_SRF,2)
            scatter3(sp(i),lsri_neutral(1,1,:),lsri_neutral(2,1,:),...
                LMSRI(plt.N1(i),j,:)./LMSRI(plt.N2(i),j,:),...
                d.s,'filled','MarkerFaceAlpha',d.MFA)
        end
        scatter3(sp(i),lsri_neutral(1,1,:),lsri_neutral(2,1,:),...
            LMSRI_neutral(plt.N1(i),1,:)./LMSRI_neutral(plt.N2(i),1,:),...
            d.s,'k','filled','MarkerFaceAlpha',d.MFA)
        xlabel('{\itl}_{MB}');
        ylabel('{\its}_{MB}');
        zlabel(sprintf('%s/%s',plt.names{plt.N1(i)},plt.names{plt.N2(i)}));
        %axis fill; grid off
        view(sp(i),[0,0]);
        axis manual
        xticks([])
        yticks([])
    end
    
    plt.auto_rotate =    0;
    plt.saveGif =        0; %save a 360 gif?
    
    if plt.auto_rotate
        spInd = [1,4,5,7,8,9];
        k=1;
        while k<360
            for j = spInd
                camorbit(sp(j),1,0,'data',[0 0 1]);
            end
            drawnow
            if plt.saveGif
                filename = [base,'\allComboSignals_rotating.gif'];% - 'signal combination'
                frame = getframe(gcf);
                im = frame2im(frame);
                [A,map] = rgb2ind(im,256);
                if k == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.005);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.005);
                end
                k=k+1;
            end
        end
    end
end

if plt.print
    save2pdf([base,'\allComboSignals.pdf'])
end

%% Random Data
% Repeat of above, but with random data 

LMSRI_rand = (randn(size(LMSRI,1),1,size(LMSRI,3))+5)*20;
for i=1:size(T_SPD,2)
    lsri_rand(1:2,:,i) = LMSToMacBoyn(LMSRI_rand(1:3,:,i),T_SSF(:,1:3)',T_lum');
    t_r(:,:,i)         = LMSToMacBoyn(LMSRI_rand([1,2,4],:,i),[T_SSF(:,1:2)';T_SSF(:,4)'],T_lum');
    t_i(:,:,i)         = LMSToMacBoyn(LMSRI_rand([1,2,5],:,i),[T_SSF(:,1:2)';T_SSF(:,5)'],T_lum');
end
lsri_rand(3,:,:) = t_r(2,:,:); clear t_r
lsri_rand(4,:,:) = t_i(2,:,:); clear t_i

%------------% Repeat of previous section

if plt.disp
    plt.names={'L','M','S','R','I'}; 
    plt.N1 = [1,NaN,NaN,1,2,NaN,1,2,3]; % Plot numbers
    plt.N2 = [2,NaN,NaN,3,3,NaN,5,5,5];
    
    %figure('units','normalized','outerposition',[0 0 1 1])    
    figure('units','normalized','outerposition',[0 0 1 1],'Renderer','opengl')
    
    for i=1:9
        if isnan(plt.N1(i))
            continue
        end
        sp(i)=subplot(3,3,i);
        hold on
        
        scatter3(sp(i),lsri_rand(1,1,:),lsri_rand(2,1,:),...
            LMSRI_rand(plt.N1(i),1,:)./LMSRI_rand(plt.N2(i),1,:),...
            d.s,'b','filled','MarkerFaceAlpha',d.MFA)
        xlabel('{\itl}_{MB}');
        ylabel('{\its}_{MB}');
        zlabel(sprintf('%s/%s',plt.names{plt.N1(i)},plt.names{plt.N2(i)}));
        %axis fill; grid off
        view(sp(i),[0,0]);
        axis manual
        xticks([])
        yticks([])
    end
    
    plt.auto_rotate =    0;
    plt.saveGif =        0; %save a 360 gif?
    
    if plt.auto_rotate
        spInd = [1,4,5,7,8,9];
        k=1;
        while k<360
            for j = spInd
                camorbit(sp(j),1,0,'data',[0 0 1]);
            end
            drawnow
            if plt.saveGif
                filename = [base,'\allComboSignals_rotating_rand.gif'];% - 'signal combination'
                frame = getframe(gcf);
                im = frame2im(frame);
                [A,map] = rgb2ind(im,256);
                if k == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.005);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.005);
                end
                k=k+1;
            end
        end
    end
end

if plt.print
    save2pdf([base,'\allComboSignals_rand.pdf'])
end



