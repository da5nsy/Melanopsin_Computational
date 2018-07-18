clear, clc, close all

% A fresh attempt

%% Load Reflectances
load sur_vrhel
refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
T_vrhel=sur_vrhel(:,refs)';
T_vrhel=T_vrhel([2,1,3,6,5,4,8,9,11,10,7],:); %change order
clear sur_vrhel refs

%% Compute Daylight SPD

D_CCT=1./linspace(1/3600,1/25000,20); %non-linear range, aiming to better reproduce observed variation
load B_cieday
T_Dspd = GenerateCIEDay(D_CCT,[B_cieday]); %these appear to be linearly upsampled from 10nm intervals (see 'cieday investigation.m' https://github.com/da5nsy/General-Purpose-Functions/blob/3ee587429e9c4f3dd52d64acd95acf82d7e05f47/cieday%20investigation.m)
T_Dspd = (T_Dspd./repmat(max(T_Dspd),81,1)); %normalise
T_Dspd = SplineSpd(S_cieday,T_Dspd,S_vrhel,2)'; % '2' flag -> Linear interp, linear extend

%% Observer

% Smith-Pokorny, for use with MacLeod Boynton diagram
load T_cones_sp %should use a different version (see 'PTB_SP.m': https://github.com/da5nsy/General-Purpose-Functions/blob/3ee587429e9c4f3dd52d64acd95acf82d7e05f47/PTB_SP.m)
load T_rods
load T_melanopsin

% Pull them all together
T_LMSRI=[(SplineCmf(S_cones_sp,T_cones_sp,S_vrhel));...
    (SplineCmf(S_rods,T_rods,S_vrhel));...
    (SplineCmf(S_melanopsin,T_melanopsin,S_vrhel))];
S_LMSRI=S_vrhel;

%% Combine

for i=1:size(T_Dspd,1)
    T_rad(:,:,i) = T_vrhel.*T_Dspd(i,:);
    LMSRI(:,:,i) = T_rad(:,:,i)*T_LMSRI';
    ls(:,:,i)    = LMSToMacBoyn(LMSRI(:,1:3,i)');
end

figure, hold on, axis equal, xlim([0 1]), ylim([0 1])
for i=1:4:size(T_Dspd,1)
    plot(ls(1,:,i),ls(2,:,i),'o-')
end

%% Factors

% Calibrating by l or s results in hitting unity in the relevent dimension
% (Perfect calibration, naff differentiation)
ls2 = ls./[ls(2,:,:);ls(2,:,:)];
figure, hold on, axis equal, %xlim([0 1]), ylim([0 1])
for i=1:4:size(T_Dspd,1)
    plot(ls2(1,:,i),ls2(2,:,i),'o-')
end

ls3 = ls./[ls(1,:,:);ls(1,:,:)];
figure, hold on, axis equal, %xlim([0 1]), ylim([0 1])
for i=1:4:size(T_Dspd,1)
    plot(ls3(1,:,i),ls3(2,:,i),'o-')
end