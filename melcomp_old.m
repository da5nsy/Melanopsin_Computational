function [MB1_minSD,MB2_minSD,melpeak,MB1_zeroSD,MB2_zeroSD,spread,MBx_m]= melcomp_old(offset)

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

%% Set-up
if ~exist('offset','var'); clear, clc, close all; end %clears everything, unless we're inside a function

% Include reflectances?
InclReflectances=1;
% Only the natural ones?
NatOnly = 1;

% Run the null condition? (Random data)
NullCondition = 0;

%% LOAD

plot_daylight=  0;
plot_obs=       0;
plot_refs=      0;

% Load Daylight Data
load('C:\Users\cege-user\Dropbox\UCL\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
granada=final; clear final
% 300 - 1100nm, 5nm interval, unlabeled
% 2600 samples
daylight=granada(17:97,:); %match obs
S_daylight=[380,5,81];
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)

if plot_daylight
    figure, hold on;
    
    plot(SToWls(S_daylight),daylight,'LineWidth',4)
    %drawnow, pause(0.3)
    
    xlabel('Wavelength (nm)')
    ylabel('Spectral Power Distribution (W m ^-^2 nm^-^1)')
    xlim([min(SToWls(S_daylight)),max(SToWls(S_daylight))])
end

% Obs data
load('T_cones_ss10')
% load('T_cones_ss2')

if plot_obs
    figure, hold on;
    for i=1:3
        plot(SToWls(S_cones_ss10),T_cones_ss10(i,:),'LineWidth',4)
        %drawnow, pause(0.3)
    end
    xlim([380 780])
    xlabel('Wavelength (nm)')
    ylabel('Normalised Spectral Sensitivity')
end

% Mel data
load('T_melanopsin')
if exist('offset','var')
    S_melanopsin(1)=S_melanopsin(1)+offset;
end
if plot_obs
    %figure,
    plot(SToWls(S_melanopsin),T_melanopsin,'LineWidth',4)
    
    legend({'L','M','S','Mel'})
end

if InclReflectances
    % Spectral Reflection Functions
    load sur_vrhel
    
    if NatOnly
        refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
        sur_vrhel=sur_vrhel(:,refs);
    end
    
    if plot_refs
        figure,
        plot(SToWls(S_vrhel),sur_vrhel,'LineWidth',4)
        if NatOnly
            legend(string(refs))
        end
    end
end

%% Initial Calculations

% Where InclReflectances is on, it fills in matrices but leaves the first
% space in each matrix for the EE ill, filled in after

if InclReflectances
    for i=1:size(sur_vrhel,2)
        
        %Interpolate to match daylight
        SFR(:,i)= interp1(SToWls(S_vrhel),sur_vrhel(:,i),SToWls(S_daylight),'linear','extrap');
        
        % %Check interpolation
        % figure, hold on
        % scatter(SToWls(S_vrhel),sur_vrhel(:,134),'r','filled')
        % scatter(SToWls(S_daylight),SFR,'g','filled')
        
        %Factor daylight by SFR
        spectra(:,:,i+1) =  repmat(SFR(:,i),1,size(daylight,2)).*daylight(:,:,1);
        
        % %Plot SFR, daylight1 and daylight 2, same graph
        % figure, hold on
        % plot(SToWls(S_daylight),SFR,'r')
        % plot(SToWls(S_daylight),daylight(:,1)*(1/(max(daylight(:,1)))),'g')
        % plot(SToWls(S_daylight),daylight2(:,1)*(1/(max(daylight2(:,1)))),'b')
        % legend({'SFR','Daylight','Daylight 2'})
        %
        % %Plot SFR, daylight1 and daylight 2, seperate graph
        % figure,
        % plot(SToWls(S_daylight),SFR,'r');ylim([0 1]);
        % figure,
        % plot(SToWls(S_daylight),daylight(:,1),'g');ylim([0 .01]);
        % figure,
        % plot(SToWls(S_daylight),daylight2(:,1),'b');ylim([0 .01]);
        
        % Calculate LMS of daylight samples
        LMS(:,:,i+1)=SplineSpd(S_cones_ss10,T_cones_ss10',S_daylight)'*spectra(:,:,i+1);
        
        % Calculate MacLeod - Boynton chromaticities of daylight samples
        MB(:,:,i+1)=LMSToMacBoyn(LMS(:,:,i+1));
        %figure,scatter(daylight_MB(1,:),daylight_MB(2,:),'k.')
        
        % Calulate M of daylight samples
        Mel(:,:,i+1)=SplineSpd(S_melanopsin,T_melanopsin',S_daylight)'*spectra(:,:,i+1);
        
    end
    
end
% Calculate LMS of daylight samples
LMS(:,:,1)=SplineSpd(S_cones_ss10,T_cones_ss10',S_daylight)'*daylight;
if NullCondition; LMS=(randn(3,500)+5)*20; end

% Calculate MacLeod - Boynton chromaticities of daylight samples
MB(:,:,1)=LMSToMacBoyn(LMS(:,:,1));
%figure,scatter(daylight_MB(1,:),daylight_MB(2,:),'k.')

% Calulate M of daylight samples
Mel(:,:,1)=SplineSpd(S_melanopsin,T_melanopsin',S_daylight)'*daylight(:,:,1);
if NullCondition; Mel=(randn(1,500)+5)*20; end

% % Calculate LMS of daylight samples
% LMS=SplineSpd(S_cones_ss10,T_cones_ss10',S_daylight)'*daylight;
%
% % Calculate MacLeod - Boynton chromaticities of daylight samples
% MB=LMSToMacBoyn(LMS);
% %figure,scatter(daylight_MB(1,:),daylight_MB(2,:),'k.')
%
% % Calulate M of daylight samples
% mel=SplineSpd(S_melanopsin,T_melanopsin',S_daylight)'*daylight;

LMSM=[LMS;Mel];

%Display Settings
dS=15;
dMEC=[.2 .2 .2];
dMFC=[.8 .8 .9];
dLW=.1;

%% What is the correlation between Mel and L/M/S?

% There is a strong correlation between Mel and all basic signals.
% This is unsurprsing as the first principal component daylight is very
% broad.

plot_corr=   0;
if plot_corr
    plotOrder={'L','M','S','Mel'};
    figure('units','normalized','outerposition',[0 0 1 1])
    
    for i=1:4
        sp(i)=subplot(2,2,i);
        scatter(...
            sp(i),...
            LMSM(i,:),...
            LMSM(4,:),...
            dS,...
            'MarkerEdgeColor',dMEC,...
            'MarkerFaceColor',dMFC,...
            'LineWidth',dLW...
            )
        xlabel(sp(i),plotOrder{i});ylabel('Mel');
        
    end
    set(subplot(2,2,4),'Color',[.8,.8,.8])
end

% Calculate r^2 values?

%% Do any signals predict MB chromaticity?

% No. They all suck at it.
% They all flatline as chromaticity changes, and then shoot up and slightly
% back on themselves in that boomerang shape.

plot_predict=   1;
if plot_predict
    plotOrder={'L','M','S','Mel'};
    
    figure('units','normalized','outerposition',[0 0 1 1])
    
    for i=1:4
        sp(i)=subplot(2,2,i);
        scatter3(...
            sp(i),...
            MB(1,:),...
            MB(2,:),...
            LMSM(i,:),...
            dS,...
            'MarkerEdgeColor',dMEC,...
            'MarkerFaceColor',dMFC,...
            'LineWidth',dLW...
            )
        
        xlim([0 1]);ylim([0 1]);%zlim([0 20])
        xlabel(sp(i),'MB1');ylabel(sp(i),'MB2');
        zlabel(sp(i),plotOrder{i});
        %view(sp(i),[70,16]);
        
    end
end
%% Does any combination of the above perform better? (Yes)

% In the following graphs, I ask whether a ratio of any of the available
% signals against any of the other available signals improves the ability
% to signal chromaticity as a one dimensional variable.

plot_comb=  0;
if plot_comb
    plotOrder={'L','M','S','Mel'};
    ScalePlot=0;
    
    figure('units','normalized','outerposition',[0 0 1 1])
    
    for i=1:16
        sp(i)=subplot(4,4,i);
        scatter3(...
            sp(i),...
            MB(1,:),...
            MB(2,:),...
            LMSM(ceil(i/4),:)./LMSM(mod(i-1,4)+1,:),...
            dS,...
            'MarkerEdgeColor',dMEC,...
            'MarkerFaceColor',dMFC,...
            'LineWidth',dLW...
            )
        % xlabel(sp(i),'MB1');ylabel(sp(i),'MB2');
        zlabel(sprintf('%s/%s',plotOrder{ceil(i/4)},plotOrder{mod(i-1,4)+1}));
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
    
    if auto_rotate
        i=1;
        while i<360
            for j=1:numel(sp)
                camorbit(sp(j),1,0,'data',[0 0 1]);
            end
            drawnow
            if saveGif
                filename = sprintf('MC_SC_%s.gif',datestr(now,'yymmddHHMMSS')); %SC - 'signal combination'
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
        LMSM(4,:)./MB(1,:),...
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
        LMSM(4,:)./MB(2,:),...
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
    plotOrder={'L','M','S','Mel'};
    
    figure('units','normalized','outerposition',[0 0 1 1]), hold on
    
    scatter3(...
        MB(1,:),...
        MB(2,:),...
        LMSM(1,:)./LMSM(4,:),...
        dS,...
        'MarkerEdgeColor',dMEC,...
        'MarkerFaceColor','r',...
        'LineWidth',dLW...
        )
    xlabel('MB1');ylabel('MB2');
    
    scatter3(...
        MB(1,:),...
        MB(2,:),...
        (LMSM(1,:)+LMSM(2,:))./LMSM(4,:),...
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

%% Reflectances

% % Plot all
% figure, hold on
% for i=1:170
%     plot(SToWls(S_vrhel),sur_vrhel(:,i),'k')
%     ylim([0 1])
%     drawnow; %pause(0.1)
% end
%
% % Define natural and non-natural
% nonNat=[46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;66;67;68;70;71;72;73;74;75;76;77;78;79;80;156;157;158;159;160;161;162;163;164;165;166;167;168;169;170;45;155];
% Nat=[15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;65;69;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;131;132;133;134;135;136;137;138;139;140;141;142;143;144;145;146;147;148;149;150;151;152;153;154];
%
% figure, hold on;
% for i=[nonNat]
%     plot(SToWls(S_vrhel),sur_vrhel(:,i),'k')
%     ylim([0 1]);
% end
% title('nonNat')
%
% figure, hold on;
% for i=[Nat]
%     plot(SToWls(S_vrhel),sur_vrhel(:,i),'k')
%     ylim([0 1]);
% end
% title('Nat')

% Picking skin colour data
% 87 - 'skin -- caucasian'
% 93 - 'skin -- African American'
% 94 - 'skin -- Asian'

% figure, hold on;
% for i=83:117
%     plot(SToWls(S_vrhel),sur_vrhel(:,i),'k','LineWidth',1)
%     ylim([0 1])
% end
% for i=[87,93,94]
%     plot(SToWls(S_vrhel),sur_vrhel(:,i),'LineWidth',4)
%     ylim([0 1])
% end

% Coloured objects - foliage or fruit
% 134	apple yellow delicious
% 137	peach skin -- ywllow
% 138	peach skin -- red
% 65	banana yellow (just turned)
% 69	ribe brown banana - [excluded as appears very smooth, possibly
% erroneous data]
% 19	Tree leaf
% 24	Bush fern-like leaf
% 140	cabbage
% 141	lettuce

%Plot all chosen SFRs
% This is redundant because similar is included in an earlier section
plot_chosen_refs=   0;
if plot_chosen_refs
    figure, hold on;
    for i=1:length(refs)
        plot(SToWls(S_vrhel),sur_vrhel(:,i),'LineWidth',4)
        
    end
    xlim([min(SToWls(S_vrhel)),max(SToWls(S_vrhel))])
    ylim([0 1])
    xlabel('Wavelength (nm)')
    ylabel('Relative Spectral Reflectance')
end

%% Basic MB plot

plot_MBbas= 1;

if plot_MBbas
    figure, hold on
    % Plot spectral locus in MB space
    
    LMS_locus=SplineSpd(S_cones_ss10,T_cones_ss10',S_daylight);
    MB_locus=LMSToMacBoyn(LMS_locus');
    fill([MB_locus(1,3:end),MB_locus(1,3)],[MB_locus(2,3:end),MB_locus(2,3)],'k','LineStyle','none','FaceAlpha','0.1')
    text_nm=string(SToWls(S_daylight))';
    text(MB_locus(1,[3:26,27,30,36,41:4:52]),MB_locus(2,[3:26,27,30,36,41:4:52]),text_nm([3:26,27,30,36,41:4:52]))
    
    
    for i = 2:size(spectra,3)
        scatter(...
            MB(1,:,i),...
            MB(2,:,i),...
            'filled')
    end
    axis equal
    xlim([0 1]);ylim([0 1]);
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
    %     LMSM(4,:,i)./(LMSM(1,:,i)+LMSM(2,:,i)),...
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


%% Hard code factors
plot_hc=    0;

if plot_hc
    fac1=   0.25;
    fac2=   -1.4;
    
    MB2=MB;
    MB2(1,:,:)=MB(1,:,:)+fac1*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:)));
    MB2(2,:,:)=MB(2,:,:)+fac2*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:)));

%     % multiplicative version
%     fac1=   -0.23;
%     fac2=   1.7;
%     
%     MB2=MB;
%     MB2(1,:,:)=MB(1,:,:).*(1-fac1*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:))));
%     MB2(2,:,:)=MB(2,:,:).*(1-fac2*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:))));
    
    figure, hold on
    for i = 2:size(spectra,3)
        scatter(...
            MB2(1,:,i),...
            MB2(2,:,i),...
            'filled')
    end
    axis equal
    xlim([0 1]);ylim([-0.5 0.5]);
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

%% Iterations

plot_it=    1;

if ~exist('offset','var') && plot_it; figure, hold on; end

MBx=MB;
MBx1std=[0;0]; %Initialise variables so that they can be added to
MBx2std=[0;0];

for S1= -10:0.05:10%-1:0.01:1.5
    MBx(1,:,:)=MB(1,:,:)+S1*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:)));
    MBx1std=[MBx1std,[S1;mean(std(MBx(1,:,2:end)))]];
end
MBx1std=MBx1std(:,2:end);

for S2= -10:0.05:10%0:0.04:2.5
    MBx(2,:,:)=MB(2,:,:)+S2*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:)));
    MBx2std=[MBx2std,[S2;mean(std(MBx(2,:,2:end)))]];
end
MBx2std=MBx2std(:,2:end);

[~,melpeakloc]=max(T_melanopsin);
melwavelength=SToWls(S_melanopsin);
melpeak=melwavelength(melpeakloc);

if plot_it
    %figure, hold on
    plot(MBx1std(1,:),MBx1std(2,:),'b')
    plot(MBx2std(1,:),MBx2std(2,:),'r')
    xlabel('weight of factor')
    ylabel('standard deviation')    
    title(sprintf('Mel peak-%d nm',melpeak))
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
%     plot(SToWls(S_cones_ss10),T_cones_ss10(i,:),'k')
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

%% What is the inter-object distance?

S_vals=[MBx1std(1,MB1_minSD_loc),MBx2std(1,MB2_minSD_loc)];

MBx(1,:,:)=MB(1,:,:)+S_vals(1)*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:)));
MBx(2,:,:)=MB(2,:,:)+S_vals(2)*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:)));

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
