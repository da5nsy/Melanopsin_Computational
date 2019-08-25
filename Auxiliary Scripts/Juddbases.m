% Load and plot Judd bases
% source:
% Judd, D.B., Macadam, D.L., Wyszecki, G., Budde, H.W., Condit, H.R., Henderson, S.T. and Simonds, J.L., 1964. Spectral Distribution of Typical Daylight as a Function of Correlated Color Temperature. Journal of the Optical Society of America, 54(8), p.1031.

%%
clear, clc, close all
d = DGdisplaydefaults;

data = xlsread('Juddbases.xlsx');
S_judd = WlsToS(data(:,1));
T_judd = data(:,2:end);

%%

figure, 
hold on
plot([min(SToWls(S_judd)),max(SToWls(S_judd))],[0,0],'k--','HandleVisibility','off')
plot(SToWls(S_judd),T_judd)
xlabel('Wavelength (nm)')

legend('Mean','$V_{1}$','$V_{2}$','$V_{3}$','$V_{4}$','interpreter','latex','Position',[0.58 0.33 0.17 0.23])

%% Granada data for comparison

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
T_granada = final; clear final
S_granana = [300,5,161];

T_granada = SplineSpd(S_granana,T_granada,S_judd);
S_granana = S_judd;

% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)

[granada_bases, score] = pca(T_granada','Centered',false);

%%
figure, 
hold on
rw = median(score); %relative weights 
plot([min(SToWls(S_judd)),max(SToWls(S_judd))],[0,0],'k--','HandleVisibility','off')
plot(SToWls(S_granana),granada_bases(:,1:5).*[0.1,rw(2:5)])

legend('Mean','$V_{1}$','$V_{2}$','$V_{3}$','$V_{4}$','interpreter','latex',...
    'Position',[0.54 0.49 0.17 0.23]);

% This seems essentially the same, though scaling is off
% Though - reading back through Judd et al's paper I note that they don't
% actually use PCA per se, "Details of the computational procedure are
% given by Simonds.16"

% Morris, R.H. and Morrissey, J.H., 1954. An Objective Method for Determination of Equivalent Neutral Densities of Color Film Images. II. Determination of Primary Equivalent Neutral Densities*. JOSA, 44(7), pp.530–534.


