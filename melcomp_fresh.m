clear, clc, close all

% A fresh attempt

% TO DO
% - check that lm work out the same when I calculate them without the PTB
%   function, just for lolz (seriously - to be careful)

%% Load Reflectances
load sur_vrhel
refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
%T_vrhel = sur_vrhel';
T_vrhel = sur_vrhel(:,refs)';
T_vrhel = T_vrhel([2,1,3,6,5,4,8,9,11,10,7],:); %change order
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

% For messing with hypothetical spectral sensitivity of melanopsin
%S_melanopsin(1)=S_melanopsin(1)+30;

% Pull them all together
T_LMSRI=[(SplineCmf(S_cones_sp,T_cones_sp,S_vrhel));...
    (SplineCmf(S_rods,T_rods,S_vrhel));...
    (SplineCmf(S_melanopsin,T_melanopsin,S_vrhel))];
S_LMSRI=S_vrhel;

%% Combine

plt_locus = 1;

for i=1:size(T_Dspd,1)
    T_rad(:,:,i) = T_vrhel.*T_Dspd(i,:);
    LMSRI(:,:,i) = T_rad(:,:,i)*T_LMSRI';
    ls(:,:,i)    = LMSToMacBoyn(LMSRI(:,1:3,i)');
    
    m(1,:,i)     = LMSRI(:,5,i)./(0.6373*LMSRI(:,1,i)+0.3924*LMSRI(:,2,i)); 
    % not called 'i' because I use that for all the loops
    % used the same scalars for luminance as are in the LMSToMAcBoyn
    %   function - which relate to the Smith-Pokorny fundamentals
end

figure, hold on, axis equal, xlim([0 1]), ylim([0 1]), zlim([0,1])
for i=1:size(T_Dspd,1)
    plot3(ls(1,:,i),ls(2,:,i),m(1,:,i),'o-')
end
xlabel('l'),ylabel('s'),zlabel('m');

view(188,46)

if plt_locus
    MB_locus=LMSToMacBoyn(T_cones_sp);
    %plot(MB_locus(1,:),MB_locus(2,:))
    fill([MB_locus(1,5:65),MB_locus(1,5)],[MB_locus(2,5:65),MB_locus(2,5)],'k','LineStyle','none','FaceAlpha','0.1')
end

%% Correction through rotation

lsm=[ls(1:2,:);m(:)']';

%rotation matrix
a=0.8036; %angle in radians, just eyeballed, and in one dimension
rm=...
    [1,0,0;...
    0,cos(a),-sin(a);...
    0,sin(a),cos(a)]; 

lsm_r=lsm*rm';

figure, hold on, axis equal, 
%xlim([0 1]), ylim([-1 1]), zlim([0,2])
plot3(lsm(:,1),lsm(:,2),lsm(:,3),'bo')
plot3(lsm_r(:,1),lsm_r(:,2),lsm_r(:,3),'ro')

legend({'Original','Rotated'},'Location','best')
xlabel('l'),ylabel('s2'),zlabel('m2'); %l stays the same

%% Correction through shift

lsm=[ls(1:2,:);m(:)']';
lsm_s = lsm; %shifted

lsm_s(:,3) = lsm(:,3)-0.27;
lsm_s(:,2) = lsm(:,2)-lsm_s(:,3);

figure, hold on, axis equal, 
plot3(lsm(:,1),lsm(:,2),lsm(:,3),'bo')
plot3(lsm_s(:,1),lsm_s(:,2),lsm_s(:,3),'ro')

legend({'Original','Shifted'},'Location','best')
xlabel('l'),ylabel('s2'),zlabel('m2');

