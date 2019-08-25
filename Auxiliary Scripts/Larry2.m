% Trying to understand how much noise there is which might perturb the
% illuminant estimate.

% The first 3 PCs of an L,M,S,I cloud under a single illuminant with natural
% surfaces can account for 99.9765 of the variance.

% This gives us an idea about how many points are away from the 'plane' but
% doesn't tell us how 'correct' that plane was...


clear, clc, clear

%% Load data

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_extended',...
    'SSF','SS10',...
    'lum','CIE_10',...
    'mel_offset',0);

%% Colorimetry

[LMSRI, ~] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

for i=1:size(T_SPD,2)
illN = i;%100;
LMSRI2 = LMSRI([1,2,3,5],:,illN);
%figure, scatter3(LMSRI(1,:,illN),LMSRI(2,:,illN),LMSRI(3,:,illN),'k.')
[~, ~, ~, ~, explained(:,i), ~] = pca(LMSRI2');
%disp(explained)
end

%%
m = mean(explained,2);
disp(sum(m(1:end-1)))