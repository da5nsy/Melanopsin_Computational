%

% Testing how close the values calculated within LMSToMacBoyn are to those
% listed in CIE 170-2:2015

clear, clc

[T_SPD, T_SRF, T_cones, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_1',...
    'SSF','SS10',...
    'mel_offset',0,...
    'lum','CIE_10');

% Compute colorimetry
[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_cones, T_lum, S_sh);

%%

T_cones = T_cones';
T_lum = T_lum';

LMS = LMSRI(1:3,1,1);

%%
factorsLM = (T_cones(1:2,:)'\T_lum');
factorS = 1/max(T_cones(3,:)./(factorsLM(1)*T_cones(1,:) + factorsLM(2)*T_cones(2,:)));
LMS = diag([factorsLM ; factorS])*LMS;

% Compute ls coordinates from LMS
n = size(LMS,2);
ls = zeros(2,n);
denom = [1 1 0]*LMS;
ls = LMS([1 3],:) ./ ([1 1]'*denom);

%%
factorsLM
factorS
