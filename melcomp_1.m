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

% Include reflectances?
InclReflectances = 1;
% Only the natural ones?
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


%% Initial Calculations

% Where InclReflectances is on, it fills in matrices but leaves the first
% space in each matrix for the EE ill, filled in after

if InclReflectances
    for i=1:size(T_SRF,2)
        %Factor daylight by SFR
        T_rad(:,:,i+1) =  repmat(T_SRF(:,i),1,size(T_SPD,2)).*T_SPD(:,:,1);
        
        % Calculate LMS of daylight samples
        LMSRI(:,:,i+1)=T_SSF'*T_rad(:,:,i+1);
        
        % Calculate MacLeod - Boynton chromaticities of daylight samples
        MB(:,:,i+1)=LMSToMacBoyn(LMSRI(1:3,:,i+1),T_SSF(:,1:3)',T_lum');
        %figure,scatter(daylight_MB(1,:),daylight_MB(2,:),'k.')
    end
end

% Calculate LMS of daylight samples
LMSRI(:,:,1) = T_SSF' * T_SPD;

% Random data
if RandomData
    LMSRI = (randn(size(T_SSF,2),500)+5)*20; 
end

% Calculate MacLeod - Boynton chromaticities of daylight samples
MB(:,:,1)=LMSToMacBoyn(LMSRI(1:3,:,1),T_SSF(:,1:3)',T_lum');

%% What is the correlation between I and L/M/S?

% There is a strong correlation between I and all basic signals.
% This is unsurprsing considering the first principal component of this
% daylight dataset.

plot_corr=   0;
if plot_corr
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

plot_predict=   0;
if plot_predict
    plotOrderNames = {'L','M','S','R','I'};
    plotOrderNums = [1,2,3,5]; 
    figure('Position',[plot_where plot_size])     
    for i=1:length(plotOrderNums)
        sp(i)=subplot(2,2,i);
        scatter3(...
            sp(i),...
            MB(1,:),...
            MB(2,:),...
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

plot_comb=  0;
if plot_comb
    plotOrderNames={'L','M','S','R','I'};
    plotOrderNums1 = [1,2,3,5,1,2,3,5,1,2,3,5,1,2,3,5]; 
    plotOrderNums2 = [1,1,1,1,2,2,2,2,3,3,3,3,5,5,5,5]; 
    ScalePlot=1;
    
    figure('units','normalized','outerposition',[0 0 1 1])
    
    for i=1:16
        sp(i)=subplot(4,4,i);
        scatter3(...
            sp(i),...
            MB(1,:),...
            MB(2,:),...
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

%% MB axes

plot_MBa=   0;
if plot_MBa
    clear sp
    
    fig=figure('units','normalized','outerposition',[0 0 1 1]);
    sp(1)=subplot(1,2,1);
    scatter3(...
        sp(1),...
        MB(1,:),...
        MB(2,:),...
        LMSRI(5,:)./MB(1,:),...
        dS,...
        'MarkerEdgeColor',dMEC,...
        'MarkerFaceColor',dMFC,...
        'LineWidth',dLW...
        )
    
    %axis fill; grid off
    zlabel('Mel/MB1');
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    view(sp(1),[0,0]);
    xlim([0,1]),ylim([0,1]),zlim([0,400])
    %sp(i).PlotBoxAspectRatioMode='manual';
    %sp(i).DataAspectRatioMode='manual';
    
    sp(2)=subplot(1,2,2);
    scatter3(...
        sp(2),...
        MB(1,:),...
        MB(2,:),...
        LMSRI(5,:)./MB(2,:),...
        dS,...
        'MarkerEdgeColor',dMEC,...
        'MarkerFaceColor',dMFC,...
        'LineWidth',dLW...
        )
    
    %axis fill; grid off
    zlabel('Mel/MB2');
    %set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    view(sp(2),[0,0]);
    xlim([0,1]),ylim([0,1]),zlim([0,400])
    %sp(i).PlotBoxAspectRatioMode='manual';
    %sp(i).DataAspectRatioMode='manual';
    
    auto_rotate=    0;
    saveGif=        0; %save a 360 gif?
    if auto_rotate
        for i = 1:35
            for j=1:numel(sp)
                camorbit(sp(j),10,0,'data',[0 0 1]);
            end
            drawnow
            if saveGif
                filename = sprintf('MC_MB_%s.gif',datestr(now,'yymmddHHMMSS'));
                frame = getframe(1);
                im{i} = frame2im(frame);
                [A,map] = rgb2ind(im{i},256);
                if i == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.005);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.005);
                end
            end
        end
    end
end

%% L vs L+M

plot_LvsLM= 0;
if plot_LvsLM
    plotOrderNames={'L','M','S','Mel'};
    
    figure('units','normalized','outerposition',[0 0 1 1]), hold on
    
    scatter3(...
        MB(1,:),...
        MB(2,:),...
        LMSRI(1,:)./LMSRI(5,:),...
        dS,...
        'MarkerEdgeColor',dMEC,...
        'MarkerFaceColor','r',...
        'LineWidth',dLW...
        )
    xlabel('MB1');ylabel('MB2');
    
    scatter3(...
        MB(1,:),...
        MB(2,:),...
        (LMSRI(1,:)+LMSRI(2,:))./LMSRI(5,:),...
        dS,...
        'MarkerEdgeColor',dMEC,...
        'MarkerFaceColor','b',...
        'LineWidth',dLW...
        )
    xlabel('MB1');ylabel('MB2');
    view([0,0]);
    
    legend('L/Mel','(L+M)/Mel')
    %xlim([0,1]),ylim([0,1]),zlim([0,4])
end

%% Basic MB plot

plot_MBbas= 0;

if plot_MBbas
    figure, hold on
    % Plot spectral locus in MB space
    
    MB_locus=LMSToMacBoyn(T_SSF(:,1:3)',T_SSF(:,1:3)',T_lum');
    fill([MB_locus(1,3:end),MB_locus(1,3)],[MB_locus(2,3:end),MB_locus(2,3)],'k','LineStyle','none','FaceAlpha','0.1')
    text_nm=string(SToWls(S_sh))';
    text(MB_locus(1,[3:26,27,30,36,41:4:52]),MB_locus(2,[3:26,27,30,36,41:4:52]),text_nm([3:26,27,30,36,41:4:52]))
    
    
    for i = 2:size(T_rad,3)
        scatter(...
            MB(1,:,i),...
            MB(2,:,i),...
            'filled')
    end
    xlim([0.4 1]);ylim([0 1]);
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    xlabel('MB1');ylabel('MB2');
    grid on
    
    % % Previous 3D model
    % figure, hold on
    % for i = 2:size(spectra,3)
    % scatter3(...
    %     MB(1,:,i),...
    %     MB(2,:,i),...
    %     LMSRI(5,:,i)./(LMSRI(1,:,i)+LMSRI(2,:,i)),...
    %     'filled')
    % end
    % axis equal
    % xlim([0 1]);ylim([0 1]);zlim([0 1]);
    % xticks(0:0.2:1)
    % yticks(0:0.2:1)
    % xlabel('MB1');ylabel('MB2');zlabel('Mel/(L+M)');
    % view(2)
    % grid on
    %
    % auto_rotate=0;
    % saveGif=0; %save a 360 gif?
    % if auto_rotate
    %     for i = 1:359
    %         camorbit(1,0,'data',[0 0 1]);
    %         drawnow
    %         if saveGif
    %             filename = sprintf('MC_M_%s.gif',datestr(now,'yymmddHHMMSS')); %M - 'model'
    %             frame = getframe(1);
    %             im{i} = frame2im(frame);
    %             [A,map] = rgb2ind(im{i},256);
    %             if i == 1
    %                 imwrite(A,map,filename1,'gif','LoopCount',Inf,'DelayTime',0.005);
    %             else
    %                 imwrite(A,map,filename1,'gif','WriteMode','append','DelayTime',0.005);
    %             end
    %         end
    %     end
    % end
    
end

%% Iterations

plot_it=    0;

if ~exist('offset','var') && plot_it; figure, hold on; end

MBx=MB;
MBx1std=[]; %Initialise variables so that they can be added to
MBx2std=[];

for S1= -10:0.05:10%-1:0.01:1.5
    MBx(1,:,:)=MB(1,:,:)+S1*(LMSRI(5,:,:)./(LMSRI(1,:,:)+LMSRI(2,:,:)));
    MBx1std=[MBx1std,[S1;mean(std(MBx(1,:,2:end)))]];
end

for S2= -10:0.05:10%0:0.04:2.5
    MBx(2,:,:)=MB(2,:,:)+S2*(LMSRI(5,:,:)./(LMSRI(1,:,:)+LMSRI(2,:,:)));
    MBx2std=[MBx2std,[S2;mean(std(MBx(2,:,2:end)))]];
end

if plot_it
    %figure, hold on
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


% %% Iterations
% if ~exist('offset','var'); figure, hold on; end
%
%
% for i=1:3
%     plot(SToWls(T_SSF),T_SSF(i,:),'k')
% end
%
% plot(SToWls(S_melanopsin),T_melanopsin,'b')
%
%
% xlabel('wavelength (nm)')
% ylabel('Spectral Sensitivity')
%
% [~,melpeakloc]=max(T_melanopsin);
% melwavelength=SToWls(S_melanopsin);
% melpeak=melwavelength(melpeakloc);
% title(sprintf('Mel peak-%d nm',melpeak))

%% Hard code factors
plot_hc=    0;

if plot_hc
    fac1=   MBx1std(1,MB1_minSD_loc); %0.25;
    fac2=   MBx2std(1,MB2_minSD_loc); %-1.4;
    
    MB2=MB;
    MB2(1,:,:)=MB(1,:,:)+fac1*(LMSRI(5,:,:)./(LMSRI(1,:,:)+LMSRI(2,:,:)));
    MB2(2,:,:)=MB(2,:,:)+fac2*(LMSRI(5,:,:)./(LMSRI(1,:,:)+LMSRI(2,:,:)));
    
    %     % multiplicative version
    %     fac1=   -0.23;
    %     fac2=   1.7;
    %
    %     MB2=MB;
    %     MB2(1,:,:)=MB(1,:,:).*(1-fac1*(LMSRI(5,:,:)./(LMSRI(1,:,:)+LMSRI(2,:,:))));
    %     MB2(2,:,:)=MB(2,:,:).*(1-fac2*(LMSRI(5,:,:)./(LMSRI(1,:,:)+LMSRI(2,:,:))));
    
    figure, hold on
    for i = 2:size(T_rad,3)
        scatter(...
            MB2(1,:,i),...
            MB2(2,:,i),...
            'filled')
    end
    xlim([0.4 1]);ylim([-0.5 0.5]);
    xticks(-1:0.2:1)
    yticks(-1:0.2:1)
    xlabel('MBx1');ylabel('MBx2')
    grid on
    
    % for i = 1:359
    %     camorbit(1,0,'data',[0 0 1]);
    %     drawnow
    %     frame = getframe(1);
    %     im{i} = frame2im(frame);
    %     [A,map] = rgb2ind(im{i},256);
    %     if i == 1
    %         imwrite(A,map,filename2,'gif','LoopCount',Inf,'DelayTime',0.005);
    %     else
    %         imwrite(A,map,filename2,'gif','WriteMode','append','DelayTime',0.005);
    %     end
    % end
    
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
