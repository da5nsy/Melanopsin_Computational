clear, clc, close all

%% Load Data
[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_1',...
    'SSF','SS10',...
    'mel_offset',0,...
    'lum','CIE_10');

%% Compute colorimetry


[~, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

lsri2 = lsri+0.1;

%%
figure,
scatter(lsri(1,:),lsri(2,:),'k.')

figure,
scatter(log(lsri(1,:)),log(lsri(2,:)),'k.')

%%
figure,
scatter(lsri2(1,:),lsri2(2,:),'k.')

figure,
scatter(log(lsri2(1,:)),log(lsri2(2,:)),'k.')
