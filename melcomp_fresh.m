clear, clc, close all

% A fresh attempt

%% Load Reflective Things

load sur_vrhel
refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
sur_vrhel=sur_vrhel(:,refs);

%% Load/Compute Daylight SPD
D_CCT=1./linspace(1/3600,1/25000,20); %non-linear range, aiming to better reproduce observed variation
load B_cieday
daylight_spd = GenerateCIEDay(D_CCT,[B_cieday]); %these appear to be linearly upsampled from 10nm intervals
daylight_spd = daylight_spd./repmat(max(daylight_spd),81,1); %normalise

%% Observer

% Smith-Pokorny, for use with MacLeod Boynton diagram
load T_cones_sp
load T_rods
load T_melanopsin

% Pull them all together
T_LMSRI=[(SplineCmf(S_cones_sp,T_cones_sp,S_vrhel));...
    (SplineCmf(S_rods,T_rods,S_vrhel));...
    (SplineCmf(S_melanopsin,T_melanopsin,S_vrhel))];
S_LMSRI=S_vrhel;

%% Convert to radiance
