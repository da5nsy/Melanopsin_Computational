%function []= messypup6(offset)

% Research questions:
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

% Answers:
%
% 1.    Melanopic signals, and signals from X,Y and Z tristimulus values
%       are all poor at predicting chromaticity. Ratios between signals are
%       all good at predicting chromaticity (apart from X:Y).

% Update notes:
%
% messypup3 - started 20171103 (cleanup and add psychtoolbox functions)
% messypup4 - swapped xyz for lms
% messypup5 - started thinking about reflectances, and line slopes
% messypup6 - added real object chromaticities

%% Set-up
if ~exist('offset','var'); clear, clc, close all; end

%Include reflectances?
InclReflectances=0;
%Only the natural ones?
NatOnly = 1;

%Run the null condition? (Random data)
NullCondition = 1;

%% LOAD

% Load Daylight Data
load('C:\Users\cege-user\Dropbox\UCL\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
granada=final; clear final
% 300 - 1100nm, 5nm interval, unlabeled
% 2600 samples
daylight=granada(17:97,1:500); %match obs
S_daylight=[380,5,81];
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hern�ndez-Andr�s, J. Romero& R.L. Lee, Jr., "Colorimetric and 
%       spectroradiometric characteristics of narrow-field-of-view 
%       clear skylight in Granada, Spain" (2001)

% Obs data
load('T_cones_ss10')
% load('T_cones_ss2')
% S_cones_ss10=S_cones_ss2;
% T_cones_ss10=T_cones_ss2;

% figure, hold on;
% for i=1:3
%     plot(SToWls(S_cones_ss10),T_cones_ss10(i,:))
%     pause(1)
% end

% Mel data
load('T_melanopsin')
if exist('offset','var')
    S_melanopsin(1)=S_melanopsin(1)+offset;
% figure,
% plot(SToWls(S_melanopsin),T_melanopsin)
end

if InclReflectances
    % Spectral Reflection Functions
    load sur_vrhel
        
    if NatOnly
        refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
        sur_vrhel=sur_vrhel(:,refs);
    end
end

%% Initial Calculations

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

% %% What is the correlation between Mel and L/M/S?
% 
% % There is a strong correlation between Mel and all basic signals.
% % This is unsurprsing as the first principal component daylight is very 
% % broad.
% 
% plotOrder={'L','M','S','Mel'};
% figure('units','normalized','outerposition',[0 0 1 1])
% 
% for i=1:4
% sp(i)=subplot(2,2,i);
% scatter(...
%     sp(i),...
%     LMSM(i,:),...
%     LMSM(4,:),...
%     dS,...
%         'MarkerEdgeColor',dMEC,...
%         'MarkerFaceColor',dMFC,...
%         'LineWidth',dLW...
%         )
% xlabel(sp(i),plotOrder{i});ylabel('Mel');
% 
% end
% set(subplot(2,2,4),'Color',[.8,.8,.8])
% 
% %% Do any signals predict MB chromaticity?
% 
% % No. They all suck at it.
% % They all flatline as chromaticity changes, and then shoot up and slightly
% % back on themselves in that boomerang shape.
% 
% plotOrder={'L','M','S','Mel'};
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% 
% for i=1:4
% sp(i)=subplot(2,2,i);
% scatter3(...
%     sp(i),...
%     MB(1,:),...
%     MB(2,:),...
%     LMSM(i,:),...
%     dS,...
%     'MarkerEdgeColor',dMEC,...
%     'MarkerFaceColor',dMFC,...
%     'LineWidth',dLW...
%     )
% 
% xlim([0 1]);ylim([0 1]);%zlim([0 20])
% xlabel(sp(i),'MB1');ylabel(sp(i),'MB2');
% zlabel(sp(i),plotOrder{i});
% %view(sp(i),[70,16]);
%  
% end
%% Does any combination of the above perform better? (Yes)

% In the following graphs, I ask whether a ratio of any of the available
% signals against any of the other available signals improves the ability
% to signal chromaticity as a one dimensional variable.

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

saveGif=0; %save a 360 gif?
i=1;
while i<360
    for j=1:numel(sp)
        camorbit(sp(j),1,0,'data',[0 0 1]);
    end
    drawnow
    if saveGif
        filename = 'signalCombination.gif';
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

% %% MB axes
% clear sp
% 
% filename = 'signalCombination.gif';
% 
% fig=figure('units','normalized','outerposition',[0 0 1 1]);
% sp(1)=subplot(1,2,1);
% scatter3(...
%     sp(1),...
%     MB(1,:),...
%     MB(2,:),...
%     LMSM(4,:)./MB(1,:),...
%     dS,...
%     'MarkerEdgeColor',dMEC,...
%     'MarkerFaceColor',dMFC,...
%     'LineWidth',dLW...
%     )
% 
% %axis fill; grid off
% zlabel('Mel/MB1');
% set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
% view(sp(1),[0,0]);
% xlim([0,1]),ylim([0,1]),zlim([0,400])
% %sp(i).PlotBoxAspectRatioMode='manual';
% %sp(i).DataAspectRatioMode='manual';
% 
% sp(2)=subplot(1,2,2);
% scatter3(...
%     sp(2),...
%     MB(1,:),...
%     MB(2,:),...
%     LMSM(4,:)./MB(2,:),...
%     dS,...
%     'MarkerEdgeColor',dMEC,...
%     'MarkerFaceColor',dMFC,...
%     'LineWidth',dLW...
%     )
% 
% %axis fill; grid off
% zlabel('Mel/MB2');
% %set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
% view(sp(2),[0,0]);
% xlim([0,1]),ylim([0,1]),zlim([0,400])
% %sp(i).PlotBoxAspectRatioMode='manual';
% %sp(i).DataAspectRatioMode='manual';
% 
% 
% for i = 1:35
%     for j=1:numel(sp)
%         camorbit(sp(j),10,0,'data',[0 0 1]);
%     end
%     drawnow    
%     frame = getframe(1);
%     im{i} = frame2im(frame);    
%     [A,map] = rgb2ind(im{i},256);
%     if i == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.005);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.005);
%     end
% end

% %% L vs L+M
% 
% plotOrder={'L','M','S','Mel'};
% 
% figure('units','normalized','outerposition',[0 0 1 1]), hold on
% 
% scatter3(...
%     MB(1,:),...
%     MB(2,:),...
%     LMSM(1,:)./LMSM(4,:),...
%     dS,...
%     'MarkerEdgeColor',dMEC,...
%     'MarkerFaceColor','r',...
%     'LineWidth',dLW...
%     )
% xlabel('MB1');ylabel('MB2');
% 
% scatter3(...
%     MB(1,:),...
%     MB(2,:),...
%     (LMSM(1,:)+LMSM(2,:))./LMSM(4,:),...
%     dS,...
%     'MarkerEdgeColor',dMEC,...
%     'MarkerFaceColor','b',...
%     'LineWidth',dLW...
%     )
% xlabel('MB1');ylabel('MB2');
% view([0,0]);
% 
% legend('L/Mel','(L+M)/Mel')
% xlim([0,1]),ylim([0,1]),zlim([0,4])

%% Reflectances

% figure, hold on
% for i=1:170
%     plot(SToWls(S_vrhel),sur_vrhel(:,i),'k')
%     ylim([0 1])
%     drawnow; %pause(0.1)
% end
 
% nonNat=[46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;66;67;68;70;71;72;73;74;75;76;77;78;79;80;156;157;158;159;160;161;162;163;164;165;166;167;168;169;170;45;155];
% Nat=[15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;65;69;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;131;132;133;134;135;136;137;138;139;140;141;142;143;144;145;146;147;148;149;150;151;152;153;154];

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
% figure, hold on;
% for i=1:length(refs)
%     plot(SToWls(S_vrhel),sur_vrhel(:,i),'LineWidth',4)
%     ylim([0 1])
% end

% %% Create scaled MB (model?)
% 
% close all
% filename1 = 'model1.gif';
% figure, hold on
% for i = 2:size(spectra,3)
% scatter3(...
%     MB(1,:,i),...
%     MB(2,:,i),...
%     LMSM(4,:,i)./(LMSM(1,:,i)+LMSM(2,:,i)),...
%     'filled')
% end
% xlim([0 1]);ylim([0 2]);zlim([0 1]);
% xlabel('MB1');ylabel('MB2');zlabel('Mel/(L+M)');
% view(3)
% grid on
% 
% for i = 1:359
%     camorbit(1,0,'data',[0 0 1]);
%     drawnow
%     frame = getframe(1);
%     im{i} = frame2im(frame);    
%     [A,map] = rgb2ind(im{i},256);
%     if i == 1
%         imwrite(A,map,filename1,'gif','LoopCount',Inf,'DelayTime',0.005);
%     else
%         imwrite(A,map,filename1,'gif','WriteMode','append','DelayTime',0.005);
%     end
% end
% 
% 
% close all
% 
% filename2 = 'model2.gif';
% MB2=MB;
% MB2(1,:,:)=MB(1,:,:)+.2*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:)));
% MB2(2,:,:)=MB(2,:,:)-1.4*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:)));
% 
% figure, hold on
% for i = 2:size(spectra,3)
% scatter3(...
%     MB2(1,:,i),...
%     MB2(2,:,i),...
%     LMSM(4,:,i)./(LMSM(1,:,i)+LMSM(2,:,i)),...
%     'filled')
% end
% xlabel('MBx1');ylabel('MBx2');zlabel('Mel/(L+M)');
% xlim([-1 1]);ylim([-1 1]);zlim([0 1]);
% view(3)
% grid on
% 
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


% %% Iterations
% if ~exist('offset','var'); figure, hold on; end
% 
% MBx=MB;
% MBx1std=[0;0];
% MBx2std=[0;0];
%  
% for S1= -10:0.05:10%-1:0.01:1.5
%     MBx(1,:,:)=MB(1,:,:)+S1*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:)));
%     MBx1std=[MBx1std,[S1;mean(std(MBx(1,:,2:end)))]];
% end
% MBx1std=MBx1std(:,2:end);
% %figure, hold on
% plot(MBx1std(1,:),MBx1std(2,:),'b')
% 
% for S2= -10:0.05:10%0:0.04:2.5
%     MBx(2,:,:)=MB(2,:,:)-S2*(LMSM(4,:,:)./(LMSM(1,:,:)+LMSM(2,:,:)));
%     MBx2std=[MBx2std,[S2;mean(std(MBx(2,:,2:end)))]];
% end
% MBx2std=MBx2std(:,2:end);
% plot(MBx2std(1,:),MBx2std(2,:),'r')
% 
% xlabel('weight of factor')
% ylabel('standard deviation')
% 
% [~,melpeakloc]=max(T_melanopsin);
% melwavelength=SToWls(S_melanopsin);
% melpeak=melwavelength(melpeakloc);
% 
% title(sprintf('Mel peak-%d nm',melpeak))

%% Iterations
if ~exist('offset','var'); figure, hold on; end


for i=1:3
    plot(SToWls(S_cones_ss10),T_cones_ss10(i,:),'k')
end

plot(SToWls(S_melanopsin),T_melanopsin,'b')


xlabel('wavelength (nm)')
ylabel('Spectral Sensitivity')

[~,melpeakloc]=max(T_melanopsin);
melwavelength=SToWls(S_melanopsin);
melpeak=melwavelength(melpeakloc);
title(sprintf('Mel peak-%d nm',melpeak))
