clc, clear, close all

%% SPDs (Spectral Power Distributions)

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)
T_SPD=final; clear final
S_SPD=[300,5,161];

%% SRFs (Spectral Reflectance Functions)

load sur_vrhel
T_SRF = sur_vrhel';
S_SRF = S_vrhel;
clear sur_vrhel S_vrhel

%% SSF (Spectral Sensitivity Functions)

% change for ss at some point in future !!!!!!!!!!1
load T_cones_sp T_cones_sp S_cones_sp
T_SSF = T_cones_sp;
S_SSF = S_cones_sp;
clear T_cones_sp S_cones_sp

load T_rods T_rods S_rods
load T_melanopsin T_melanopsin S_melanopsin
T_mel = SplineCmf(S_melanopsin,T_melanopsin, S_melanopsin - [10, 0, 0],1); %Increasing the range of this function in case it ends up limiting the range of S_sh, and shorten variable names
S_mel = S_melanopsin - [10, 0, 0]; clear S_melanopsin T_melanopsin

%% Interpolate functions to match range and interval

%reduce all data down to the common range/interval
S_sh = [max([S_SPD(1),S_SRF(1),S_SSF(1)]),max([S_SPD(2),S_SRF(2),S_SSF(2)]),min([S_SPD(3),S_SRF(3),S_SSF(3)])]; %S_shared: work out what the lowest common denominator for the range/interval of the data is

T_SPD = SplineSpd(S_SPD,T_SPD,S_sh,1); % extend == 1: Cubic spline, extends with last value in that direction
T_SRF = SplineSrf(S_SRF,T_SRF',S_sh,1); 
T_SSF  = SplineCmf(S_SSF,T_SSF,S_sh,1)';
T_rods = SplineCmf(S_rods,T_rods,S_sh)'; %extend? !!!!!!!!!!
T_mel  = SplineCmf(S_mel,T_mel,S_sh)'; %extend? !!!!!!!!!!
[S_SPD, S_SRF, S_SSF, S_rods, S_mel] = deal(S_sh);

% combine sensitivity vectors
T_LMSRI=[T_SSF,T_rods,T_mel];
S_LMSRI=S_sh;

%% Calculate colour signals

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_SRF.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_LMSRI'*T_rad(:,:,i);
    lsri(1:2,:,i) = LMSToMacBoyn(LMSRI(1:3,:,i)); %change this to be more direct !!!!!!!!!!!!!!   
    lsri(3,:,i)   = LMSRI(4,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); %change this to be more direct !!!!!!!!!!!!!!   
    lsri(4,:,i)   = LMSRI(5,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); %change this to be more direct !!!!!!!!!!!!!!   
end

%% PCA of illums

[pc.COEFF, pc.SCORE, pc.LATENT, pc.TSQUARED, pc.EXPLAINED] = pca(T_SPD'); %correct/weight at edges of spectrum?

%% Calc correlation between 

cs = cat(1,LMSRI,lsri,lsri(2,:,:)./lsri(4,:,:),lsri(2,:,:)./lsri(3,:,:)); %colour signals

for j=1:6
    for i=1:size(cs,1)
        c(i,j) = corr(squeeze(mean(cs(i,:,:),2)),pc.SCORE(:,j));
    end
end

c_norm = c;
for i=1:j %leftover, be careful
    c_norm(:,i) = c_norm(:,i) - min(c_norm(:,i));
    c_norm(:,i) = c_norm(:,i)/max(c_norm(:,i));
end

corr(squeeze(mean(lsri(2,:,:),2))./squeeze(mean(lsri(4,:,:),2)),pc.SCORE(:,2)); %s/i
surf(c_norm)





