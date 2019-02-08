function [MB1_minSD,MB2_minSD,MB1_zeroSD,MB2_zeroSD,spread,MBx_m]= melcomp_1(offset)

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
% Debug melomp_caller (try to reproduce original performance / work out why
%   it has changed.
% Make label axes italics where appropriate
% Add save commands
% Set limits for MB axes in line with other scripts

%% Pre-flight

try
    nargin;
catch
    clear, clc, close all;
    offset = 0;
end

% Only the natural reflectances?
NatOnly = 1;

% Run the null condition? (Random data)
RandomData = 0;

% Display Settings
dS=15;
dMEC=[.2 .2 .2];
dMFC=[.8 .8 .9];
dLW=.1;

plot_where = [500,200];
plot_size  = [800,400];

set(0,'defaultAxesFontName', 'Courier')

base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\melcomp_1';
disp_figures = 1;
print_figures = 0;

min_l_scale = 0.62;
max_l_scale = 0.82;
max_s_scale = 0.04;

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

[LMSRI_neutral, lsri_neutral] = melcomp_colorimetry(T_SPD, ones(S_sh(3),1), T_SSF, T_lum, S_sh);

if RandomData
    LMSRI = rand(size(LMSRI))*20;
    lsri  = rand(size(lsri)) *20;
end


%% What is the correlation between I and L/M/S?

% There is a strong correlation between I and all basic signals.
% This is unsurprsing considering the first principal component of this
% daylight dataset.

plot_corr = 1;
if plot_corr && disp_figures
    plotOrderNames = {'L','M','S','R','I'};
    plotOrderNums = [1,2,3,5];
    figure('Position',[plot_where plot_size])
    for i=1:length(plotOrderNums)
        sp(i)=subplot(2,2,i);
        scatter(...
            sp(i),...
            LMSRI(plotOrderNums(i),:),...
            LMSRI(5,:),...
            dS,...
            'MarkerEdgeColor',dMEC,...
            'MarkerFaceColor',dMFC,...
            'LineWidth',dLW...
            )
        
        xlabel(sp(i),['{\it',plotOrderNames{plotOrderNums(i)},'}']);
        ylabel('{\itI}');
        xticks([min(xlim),max(xlim)])
        yticks([min(ylim),max(ylim)])
    end
    set(subplot(2,2,4),'Color',[.8,.8,.8])
end

if print_figures
    save2pdf([base,'\correlation.pdf'])
end

%% Do any signals predict MB chromaticity?

% No. They all suck at it.
% They all flatline as chromaticity changes, and then shoot up and slightly
% back on themselves in that boomerang shape.

plot_predict = 1;
if plot_predict && disp_figures
    plotOrderNames = {'L','M','S','R','I'};
    plotOrderNums = [1,2,3,5];
    figure('Position',[plot_where plot_size])
    for i=1:length(plotOrderNums)
        sp(i)=subplot(2,2,i);
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
        
        xlim([min_l_scale max_l_scale]);
        ylim([0 max_s_scale]);
        xlabel('{\itl}_{MB}');
        ylabel('{\its}_{MB}');
        zlabel(sp(i),['{\it',plotOrderNames{plotOrderNums(i)},'}']);
        %view(sp(i),[70,16]);
        xticks([min(xlim),max(xlim)])
        yticks([min(ylim),max(ylim)])
        zticks([min(zlim),max(zlim)])
    end
end

if print_figures
    save2pdf([base,'\predict.pdf'])
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
    
    auto_rotate=    1;
    saveGif=        0; %save a 360 gif? %Not currently working !!!!!!!!!!!!!!!!!!!!!!!1
    
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

%% Basic MB plot

plot_MBbas = 1;

if plot_MBbas && disp_figures
    figure('Position',[plot_where plot_size]), hold on
    %MB_locus=LMSToMacBoyn(T_SSF(:,1:3)',T_SSF(:,1:3)',T_lum');
    %fill([MB_locus(1,3:end),MB_locus(1,3)],[MB_locus(2,3:end),MB_locus(2,3)],'k','LineStyle','none','FaceAlpha','0.1')
    %text_nm=string(SToWls(S_sh))';
    %text(MB_locus(1,[3:26,27,30,36,41:4:52]),MB_locus(2,[3:26,27,30,36,41:4:52]),text_nm([3:26,27,30,36,41:4:52]))
    
    scatter(lsri(1,:),lsri(2,:),'filled')
    
    xlim([min_l_scale max_l_scale]);ylim([0 max_s_scale]);
    xticks([min(xlim),max(xlim)])
    yticks([min(ylim),max(ylim)])
    xlabel('{\itl}_{MB}');
    ylabel('{\its}_{MB}');
end

%% Iterations

plot_it = 1;

if plot_it && disp_figures
    figure('Position',[plot_where plot_size]), hold on
end

MBx=lsri;
MBx1std=[]; %Initialise variables so that they can be added to
MBx2std=[];

for S1= -3:0.01:3
    MBx(1,:,:)=lsri(1,:,:)+S1*lsri(4,:,:);
    MBx1std=[MBx1std,[S1;mean(std(MBx(1,:,:)))]];
end

for S2= -3:0.01:3
    MBx(2,:,:)=lsri(2,:,:)+S2*lsri(4,:,:);
    MBx2std=[MBx2std,[S2;mean(std(MBx(2,:,:)))]];
end

if plot_it
    plot(MBx1std(1,:),MBx1std(2,:),'b')
    plot(MBx2std(1,:),MBx2std(2,:),'r')
    xlabel('weight of factor')
    ylabel('standard deviation')
    title(sprintf('Mel peak-%d nm',488+offset))
    legend({'MBx1','MBx2'})
end

[MB1_minSD, MB1_minSD_loc] = min(MBx1std(2,:));
[MB2_minSD, MB2_minSD_loc] = min(MBx2std(2,:));
MB1_zeroSD = MBx1std(2,MBx1std(1,:)==0);
MB2_zeroSD = MBx2std(2,MBx2std(1,:)==0);

%% Hard code factors
plot_hc = 1;

if plot_hc && disp_figures
    fac1=   MBx1std(1,MB1_minSD_loc);
    fac2=   MBx2std(1,MB2_minSD_loc);
    
    MB_star=lsri(1:2,:,:);
    MB_star(1,:,:)=lsri(1,:,:) + fac1 * lsri(4,:,:);
    MB_star(2,:,:)=lsri(2,:,:) + fac2 * lsri(4,:,:);
    
    figure('Position',[plot_where plot_size]), hold on
    scatter(MB_star(1,:),MB_star(2,:),'filled')   
    
    xticks([])
    yticks([])
    xlabel('MBx1');ylabel('MBx2')
    grid on
    
end

%% What is the inter-object distance?

S_vals=[MBx1std(1,MB1_minSD_loc),MBx2std(1,MB2_minSD_loc)];

MBx(1,:,:)=MB(1,:,:)+S_vals(1)*(LMSRI(5,:,:)./(LMSRI(1,:,:)+LMSRI(2,:,:)));
MBx(2,:,:)=MB(2,:,:)+S_vals(2)*(LMSRI(5,:,:)./(LMSRI(1,:,:)+LMSRI(2,:,:)));

% for i=2:12
%     scatter(squeeze(MBx(1,:,i)),squeeze(MBx(2,:,i)));
%     xlim([-10 10])
%     ylim([-10 10])
%     hold on
%     axis equal
% end
% hold off

MBx_m=squeeze(mean(MBx(:,:,2:end),2));
% scatter(MBx_m(1,2:end),MBx_m(2,2:end),'k','filled');
% axis equal
% xlim([0 1])
% ylim([-1 2])

spread=[std(MBx_m(1,:)),std(MBx_m(2,:))];
%disp(melpeak)
