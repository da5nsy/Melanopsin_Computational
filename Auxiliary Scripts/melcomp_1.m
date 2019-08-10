function [MB1_minSD,MB2_minSD,MB1_baselineSD,MB2_baslineSD,spread,MB_star_m]= melcomp_1(offset,varargin)

%% Research questions:
%
% 1.	Considering only daylight spectra (excluding reflective surfaces),
%       can a melanopic signal predict the chromaticity of daylight more
%       precisely than signals provided by other retinal cell populations?
%
% 2.	Now considering also object reflectances, does a melanopic signal
%       provide a means of calculating a sign and weight of shift required
%       to counteract the chromatic shift induced upon objects by a change
%       in daylight conditions?
%
% 3.	Are there objects, or luminance levels, or daylight chromaticities
%       for which the melanopic signal is particularly effective or
%       ineffective at performing the above task

% TO-DO
% Fix gif saving issue
% Would be nice to decrease file/loading size of full-page graphics

%% Pre-flight

try
    nargin;
catch
    clear, clc, close all;
    offset = 0;
    varargin = {};
end

default_plt_appf_overide = 0;
default_plt_it_overide = 0;
p = inputParser;
addParameter(p,'plt_appf_overide',default_plt_appf_overide);
addParameter(p,'plt_it_overide',default_plt_it_overide);
parse(p,varargin{:});

% Only the natural reflectances?
NatOnly = 1;

% Run the Random data condition?
RandomData = 0;

% Display Settings
dS=10;
dMEC=[.2 .2 .2];
dMFC=[.8 .8 .9];
dLW=.1;

set(groot,'defaultfigureposition',[100 100 500 400]); 
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultAxesFontName', 'Courier');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultFigureRenderer', 'painters') %renders pdfs as vectors
set(groot,'defaultfigurecolor','white')

base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\melcomp_1';
disp_figures  = 1;
print_figures = 0;

min_l_scale = 0.62;
max_l_scale = 0.82;
max_s_scale = 0.05;

%% LOAD

if NatOnly
    refs = 'Vrhel_nat_1';
else
    refs = 'Vrhel_full';
end

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada',...
    'SRF',refs,...
    'SSF','SS10',...
    'lum','CIE_10',...
    'mel_offset',offset);


%% Colorimetry

[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

% generate colorimetry for the light source itself
[LMSRI_neutral, lsri_neutral] = melcomp_colorimetry(T_SPD, ones(S_sh(3),1), T_SSF, T_lum, S_sh);

%generate random data to examine inherent relationships due to col calcs
LMSRI_rand = (randn(size(LMSRI,1),1,size(LMSRI,3))+5)*20;
for i=1:size(T_SPD,2)
    lsri_rand(1:2,:,i) = LMSToMacBoyn(LMSRI_rand(1:3,:,i),T_SSF(:,1:3)',T_lum');
    t_r(:,:,i)         = LMSToMacBoyn(LMSRI_rand([1,2,4],:,i),[T_SSF(:,1:2)';T_SSF(:,4)'],T_lum');
    t_i(:,:,i)         = LMSToMacBoyn(LMSRI_rand([1,2,5],:,i),[T_SSF(:,1:2)';T_SSF(:,5)'],T_lum');
end
lsri_rand(3,:,:) = t_r(2,:,:); clear t_r
lsri_rand(4,:,:) = t_i(2,:,:); clear t_i


%% What is the correlation between I and L/M/S?

% There is a strong correlation between I and all basic signals.
% This is unsurprsing considering the first principal component of this
% daylight dataset.

plot_corr = 1;
if plot_corr && disp_figures
    plotOrderNames = {'L','M','S','R','I'};
    plotOrderNums = [1,2,3];
    figure, hold on
    for i=1:length(plotOrderNums)
        sp(i)=subplot(1,3,i); hold on
        scatter(...
            sp(i),...
            LMSRI(plotOrderNums(i),:),...
            LMSRI(5,:),...
            dS,...
            'MarkerEdgeColor',dMEC,...
            'MarkerFaceColor',dMFC,...
            'LineWidth',dLW...
            )
        
        scatter(...
            sp(i),...
            LMSRI_neutral(plotOrderNums(i),:),...
            LMSRI_neutral(5,:),...
            dS,...
            'MarkerEdgeColor',dMEC,...
            'MarkerFaceColor','r',...
            'LineWidth',dLW...
            )
        
        xlabel(sp(i),['{\it',plotOrderNames{plotOrderNums(i)},'}']);
        xticks([min(xlim),max(xlim)])
        yticks([])
    end
    ylabel(sp(1),'{\itI}');
    yticks(sp(1),[min(ylim),max(ylim)])
end

if print_figures
    save2pdf([base,'\correlationBetweenLevel1Sigs.pdf'])
end

%% Do any signals predict MB chromaticity?

% No. They all suck at it.
% They all flatline as chromaticity changes, and then shoot up and slightly
% back on themselves in that boomerang shape.

plot_predict = 1;
if plot_predict && disp_figures
    plotOrderNames = {'L','M','S','R','I'};
    plotOrderNums = [1,2,3,5];
    figure
    for i=1:length(plotOrderNums)
        sp(i)=subplot(2,2,i);
        hold on
        scatter3(...
            sp(i),...
            lsri(1,:),...
            lsri(2,:),...
            LMSRI(plotOrderNums(i),:),...
            dS,...
            'MarkerEdgeColor',dMEC,...
            'MarkerFaceColor',dMFC,...
            'LineWidth',dLW...
            )
        scatter3(...
            sp(i),...
            lsri_neutral(1,:),...
            lsri_neutral(2,:),...
            LMSRI_neutral(plotOrderNums(i),:),...
            dS,...
            'MarkerEdgeColor',dMEC,...
            'MarkerFaceColor','r',...
            'LineWidth',dLW...
            )        
        view(sp(i),[-64,11]);
        xlim([min_l_scale max_l_scale]);
        ylim([0 0.05]);
        zlim([0 ceil(max(zlim)/5)*5]); %round up to next 10
        zlabel(sp(i),['{\it',plotOrderNames{plotOrderNums(i)},'}']);
        if i==1
            xticks([min(xlim),max(xlim)])
            yticks([min(ylim),max(ylim)])
            zticks([min(zlim),max(zlim)])
            xlabel('{\itl}_{MB}');
            ylabel('{\its}_{MB}');            
        else
            xticks([])
            yticks([])
            zticks([min(zlim),max(zlim)])
        end
    end
end

if print_figures
    save2pdf([base,'\level1sigspredictingColorimetry.pdf'])
end


%% Does any combination of the above perform better? (Yes)

% In the following graphs, I ask whether a ratio of any of the available
% signals against any of the other available signals improves the ability
% to signal chromaticity as a one dimensional variable.

plot_comb = 1;
if plot_comb && disp_figures
    plotOrderNames={'L','M','S','R','I'};
    plotOrderNums1 = [1,2,3,5,1,2,3,5,1,2,3,5,1,2,3,5];
    plotOrderNums2 = [1,1,1,1,2,2,2,2,3,3,3,3,5,5,5,5];
    ScalePlot=1;
    
    figure('units','normalized','outerposition',[0 0 1 1])
    
    for i=1:16
        sp(i)=subplot(4,4,i);
        hold on
        if RandomData
            scatter3(...
                sp(i),...
                lsri_rand(1,:),...
                lsri_rand(2,:),...
                LMSRI_rand(plotOrderNums1(i),:)./LMSRI_rand(plotOrderNums2(i),:),...
                dS,...
                'MarkerEdgeColor',dMEC,...
                'MarkerFaceColor','b',...
                'LineWidth',dLW...
                )
        else
            scatter3(...
                sp(i),...
                lsri(1,:),...
                lsri(2,:),...
                LMSRI(plotOrderNums1(i),:)./LMSRI(plotOrderNums2(i),:),...
                dS,...
                'MarkerEdgeColor',dMEC,...
                'MarkerFaceColor',dMFC,...
                'LineWidth',dLW...
                )
            scatter3(...
                sp(i),...
                lsri_neutral(1,:),...
                lsri_neutral(2,:),...
                LMSRI_neutral(plotOrderNums1(i),:)./LMSRI_neutral(plotOrderNums2(i),:),...
                dS,...
                'MarkerEdgeColor',dMEC,...
                'MarkerFaceColor','r',...
                'LineWidth',dLW...
                )
        end
        % xlabel(sp(i),'MB1');ylabel(sp(i),'MB2');
        zlabel(sprintf('%s/%s',plotOrderNames{plotOrderNums1(i)},plotOrderNames{plotOrderNums2(i)}));
        %axis fill; grid off
        view(sp(i),[0,0]);
        if ~ScalePlot
            xlim([0,1]),ylim([0,1]),zlim([0,4])
        end
        %     set(gca,'XTickLabel',[]);
        %     set(gca,'YTickLabel',[]);
    end
    
    set(sp(1),'Color',[.8,.8,.8])
    set(sp(6),'Color',[.8,.8,.8])
    set(sp(11),'Color',[.8,.8,.8])
    set(sp(16),'Color',[.8,.8,.8])
    %
    
    auto_rotate=    0;
    saveGif=        0; %save a 360 gif? 
    %Not currently working !!!!!!!!!!!!!!!!!!!!!!!1
    
    if auto_rotate
        i=1;
        while i<360
            for j=1:numel(sp)
                camorbit(sp(j),1,0,'data',[0 0 1]);
            end
            drawnow
            if saveGif
                filename = [base,'\' sprintf('MC_SC_%s.gif',datestr(now,'yymmddHHMMSS'))]; %SC - 'signal combination'
                frame = getframe(1);
                im{i} = frame2im(frame);
                [A,map] = rgb2ind(im{i},256);
                if i == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.005);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.005);
                end
                i=i+1;
            end
        end
    end
end

if print_figures
    if RandomData
        save2pdf([base,'\allComboSignals_rand.pdf'])
    else
        save2pdf([base,'\allComboSignals.pdf'])
    end
end

%% Basic MB plot

plot_MBbas = 1;

if plot_MBbas && disp_figures
    figure, hold on
    %MB_locus=LMSToMacBoyn(T_SSF(:,1:3)',T_SSF(:,1:3)',T_lum');
    %fill([MB_locus(1,3:end),MB_locus(1,3)],[MB_locus(2,3:end),MB_locus(2,3)],'k','LineStyle','none','FaceAlpha','0.1')
    %text_nm=string(SToWls(S_sh))';
    %text(MB_locus(1,[3:26,27,30,36,41:4:52]),MB_locus(2,[3:26,27,30,36,41:4:52]),text_nm([3:26,27,30,36,41:4:52]))
    
    dMECalt = hsv(size(T_SRF,2));
    
    for i=1:size(T_SRF,2)
        scatter(lsri(1,i,:),lsri(2,i,:),dS,...
            'MarkerEdgeColor',dMECalt(i,:),...
            'MarkerFaceColor',dMFC,...
            'LineWidth',dLW...
            )
    end
    
    scatter(lsri_neutral(1,:),lsri_neutral(2,:),dS,...
        'MarkerEdgeColor',dMEC,...
        'MarkerFaceColor','r',...
        'LineWidth',dLW...
        )
    
    xlim([min_l_scale max_l_scale]);
    ylim([0 max_s_scale]);
    xticks([min(xlim),max(xlim)])
    yticks([min(ylim),max(ylim)])
    xlabel('{\itl}_{MB}');
    ylabel('{\its}_{MB}');
end

if print_figures
    save2pdf([base,'\BasicMB.pdf'])
end


%% Iterations

plot_it = 1;

if plot_it && disp_figures
    figure, hold on
end

MBx=lsri;
MBx1std=[]; %Initialise variables so that they can be added to
MBx2std=[];

for S1= -2:0.01:2
    MBx(1,:,:)=lsri(1,:,:)+S1*lsri(4,:,:);
    MBx1std=[MBx1std,[S1;mean(std(MBx(1,:,:)))]];
end

for S2= -2:0.01:2
    MBx(2,:,:)=lsri(2,:,:)+S2*lsri(4,:,:);
    MBx2std=[MBx2std,[S2;mean(std(MBx(2,:,:)))]];
end

if or(plot_it && disp_figures,p.Results.plt_it_overide)
    plot3(MBx1std(1,:),MBx1std(2,:),ones(length(MBx1std),1)*offset,'r')
    plot3(MBx2std(1,:),MBx2std(2,:),ones(length(MBx1std),1)*offset,'b')
    xlabel('weight of factor')
    ylabel('standard deviation')
    legend({'MBx1','MBx2'})
    if p.Results.plt_it_overide
        title(sprintf('Mel peak-%d nm',488+offset))
    end
end

[MB1_minSD, MB1_minSD_loc] = min(MBx1std(2,:));
[MB2_minSD, MB2_minSD_loc] = min(MBx2std(2,:));
MB1_baselineSD = MBx1std(2,MBx1std(1,:)==0);
MB2_baslineSD = MBx2std(2,MBx2std(1,:)==0);

if print_figures
    save2pdf([base,'\minimiseSD.pdf'])
end


%% Apply factors
plot_appf = 1;

fac1=   MBx1std(1,MB1_minSD_loc);
fac2=   MBx2std(1,MB2_minSD_loc);

MB_star=lsri(1:2,:,:);
MB_star(1,:,:)=lsri(1,:,:) + fac1 * lsri(4,:,:);
MB_star(2,:,:)=lsri(2,:,:) + fac2 * lsri(4,:,:);

if or(plot_appf && disp_figures,p.Results.plt_appf_overide)
    
    figure, hold on
    for i=1:size(T_SRF,2)
        scatter(MB_star(1,i,:),MB_star(2,i,:),dS,...
            'MarkerEdgeColor',dMECalt(i,:),...
            'MarkerFaceColor',dMFC,...
            'LineWidth',dLW...
            )
    end
    
    cleanTicks
    xlabel('{\itl}_{MB} + {\itk_1i}_{MB}');
    ylabel('{\its}_{MB} + {\itk_2i}_{MB}');
    
    grid on
end

if print_figures
    save2pdf([base,'\correctedChromaticities.pdf'])
end

%% What is the inter-object distance?

MB_star_m=squeeze(mean(MB_star,3));

spread_l = pdist(MB_star_m(1,:)');
spread_s = pdist(MB_star_m(2,:)');

spread = [mean(spread_l), mean(spread_s)];


