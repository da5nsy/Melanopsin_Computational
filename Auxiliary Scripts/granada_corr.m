% I originally thought that there was a bug in the Granada data, because of
% the lack of correlation at certain values, but it seems to be explained by
% the mean (or sd) of the data.

% There's a weird thing happening at the lowest 5 values.
% It would be nice to have the axis tick labels being the wavelength rather
% than the index number.

clear, clc, close all

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)
T_SPD=final; clear final
S_SPD=[300,5,161];

%trying to get rid of the weird first 5 values that show 0 correlation,
%even within themselves (!?) :

%T_SPD = T_SPD - min(T_SPD(:);
%or
%T_SPD(T_SPD == 0) = min(T_SPD(:));

% I think it is because they have 0 values and that throws off the
% correlation function for some reason


figure, hold on
cr = real(corr(log2(T_SPD')));
imagesc(cr)

axis image
colorbar

%figure,
plot(mean(T_SPD,2)*200,'k','LineWidth',2)

%%
% Looks like there are only zeros in those first 5 rows of data

figure,
t=T_SPD == 0; 
imshow(t(1:10,:))
axis auto