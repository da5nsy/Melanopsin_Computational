clear, clc, close all

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)
T_SPD=final; clear final
S_SPD=[300,5,161];

%%

S_SPD_e = SToWls(S_SPD);
%T_SPD_e = T_SPD(21:81,:); %400-700nm
T_SPD_e = T_SPD(17:97,:); %380-780nm

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(T_SPD_e');



%%

load B_cieday.mat

figure, hold on

plot(SToWls(S_cieday),B_cieday(:,1:3)./max(B_cieday(:,1:3)),'r:','LineWidth',3)
plot(S_SPD_e(17:97),COEFF(:,1:5))
legend


