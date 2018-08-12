clc, clear, close all

load sur_vrhel
refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141]; %natural ones
T_refs_all = sur_vrhel';
T_refs_nat = sur_vrhel(:,refs)';
S_refs = S_vrhel; clear S_vrhel

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)

% T_SPD=final; clear final
% S_SPD=[300,5,161];

%match to cones instead
T_SPD = final(17:97,:); clear final;
S_SPD = [380,5,81];

% figure, hold on;
% plot(SToWls(S_SPD),T_SPD)
% %drawnow, pause(0.3)
% xlabel('Wavelength (nm)')
% ylabel('Spectral Power Distribution (W m ^-^2 nm^-^1)')
% xlim([S_SPD(1),max(SToWls(S_SPD))]), ylim([0 max(T_SPD(:))]);

%normalise T_SPD for luminance
load T_cones_sp

SPD_L = (0.6373*T_cones_sp(1,:)*T_SPD)+(0.3924*T_cones_sp(2,:)*T_SPD);
T_SPD_n = T_SPD./SPD_L;

figure, hold on;
plot(SToWls(S_SPD),T_SPD_n)


%% Testing STD of refs

% Natural, wavelength
figure, hold on
plot(SToWls(S_refs),T_refs_nat)
plot(SToWls(S_refs),std(T_refs_nat),'ko','DisplayName','STD')

% Natural, frequency
figure, hold on
plot(1./SToWls(S_refs),T_refs_nat)
plot(1./SToWls(S_refs),std(T_refs_nat),'ko','DisplayName','STD')

% % All, wavelength
% figure, hold on
% plot(SToWls(S_refs),T_refs_all)
% plot(SToWls(S_refs),std(T_refs_all),'ko','DisplayName','STD')
% 
% % All, frequency
% figure, hold on
% plot(1./SToWls(S_refs),T_refs_all)
% plot(1./SToWls(S_refs),std(T_refs_all),'ko','DisplayName','STD')

% Natural, wavelength
figure, hold on
plot(SToWls(S_SPD),T_SPD_n)
plot(SToWls(S_SPD),std(T_SPD_n'),'ko','DisplayName','STD')

% Natural, frequency
figure, hold on
plot(1./SToWls(S_SPD),T_SPD_n)
plot(1./SToWls(S_SPD),std(T_SPD_n'),'ko','DisplayName','STD')

%%
load T_cones_ss2.mat
load T_melanopsin.mat

%T_refs_i = SplineSrf(S_refs,T_refs_nat',S_cones_ss2)';
T_refs_i = SplineSrf(S_refs,T_refs_all',S_cones_ss2)';

%c = corr(T_refs_i);
c = corr(T_SPD_n');
c(isnan(c))=1;

figure, hold on
imagesc(c)
axis image

%%
% plot locations of peak sensitvities
[~, maxloc]        = max(T_cones_ss2,[],2);
[~, maxloc(end+1)] = max(T_melanopsin,[],2);
for i=1:4
    plot([maxloc(i),maxloc(i)],[min(ylim),max(ylim)],'k')
end
colorbar

fill([0, max(xlim), 0, 0],[0, max(ylim), max(ylim), 0],[1,1,1],'LineStyle','none')

%xlim([1 341]);ylim([1 341]);

% be careful with the below, if fairly bluntly overwrites the scales and so
% it is possible that a bug could creep in here
set(gca,'XTickLabel',390:50:790)
set(gca,'YTickLabel',390:50:790)


