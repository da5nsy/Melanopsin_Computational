clear, clc, close all

% A fresh attempt

% TO DO
% - check that lm work out the same when I calculate them without the PTB
%   function, just for lolz (seriously - to be careful)
% - convert loop names so that I can use i for melanopsin variable name (OR
%   never compute it alone, bundle it with other variables)
% - what about log signals?
% - Can the 'correction through shift' be done by division rather than subtraction?

%% Load Reflectances
load sur_vrhel
refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
%T_vrhel = sur_vrhel';
T_vrhel = sur_vrhel(:,refs)';
T_vrhel = T_vrhel([2,1,3,6,5,4,8,9,11,10,7],:); %change order for clearer plotting visualisation
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

plt_fig     = 0;
plt_locus   = 0;

for i=1:size(T_Dspd,1)
    T_rad(:,:,i) = T_vrhel.*T_Dspd(i,:);
    LMSRI(:,:,i) = T_rad(:,:,i)*T_LMSRI';
    lsri(1:2,:,i)  = LMSToMacBoyn(LMSRI(:,1:3,i)');    
    lsri(3,:,i)    = LMSRI(:,4,i)./(0.6373*LMSRI(:,1,i)+0.3924*LMSRI(:,2,i)); 
    lsri(4,:,i)    = LMSRI(:,5,i)./(0.6373*LMSRI(:,1,i)+0.3924*LMSRI(:,2,i)); 
    % used the same scalars for luminance as are in the LMSToMAcBoyn
    %   function - which relate to the Smith-Pokorny fundamentals
end

% Compute colorimetry (just for display)
load T_xyz1931.mat
T_xyz1931=SplineCmf(S_xyz1931,T_xyz1931,S_vrhel);
for i=[11, 1:size(T_Dspd,1)] %starts with 11 (6551K, arbitrary), so that a fixed white is already calculated in time for the first 'plt_RGB' line
    plt_whiteXYZ(:,i) = T_Dspd(i,:) * T_xyz1931';
    plt_XYZ(:,:,i)    = T_rad(:,:,i) * T_xyz1931';
    plt_Lab(:,:,i)    = XYZToLab(squeeze(plt_XYZ(:,:,i))',plt_whiteXYZ(:,i));    
    plt_RGB(:,:,i)    = XYZToSRGBPrimary(LabToXYZ(plt_Lab(:,:,i),plt_whiteXYZ(:,11))); %Using fixed, arbitrary (mid-range), white.
    plt_RGB(:,:,i)    = plt_RGB(:,:,i)/max(max(plt_RGB(:,:,i)));
end

if plt_fig
    figure, hold on, axis equal, xlim([0 1]), ylim([0 1]), zlim([0,1])
    for i=1:size(T_Dspd,1)
        plot3(lsri(1,:,i),lsri(2,:,i),lsri(4,:,i),'k')        
        scatter3(lsri(1,:,i),lsri(2,:,i),lsri(4,:,i),[],plt_RGB(:,:,i)','v','filled')
    end
    xlabel('l'),ylabel('s'),zlabel('m');
    
    view(188,46)
    
    if plt_locus
        MB_locus=LMSToMacBoyn(T_cones_sp);
        %plot(MB_locus(1,:),MB_locus(2,:))
        fill([MB_locus(1,5:65),MB_locus(1,5)],[MB_locus(2,5:65),MB_locus(2,5)],'k','LineStyle','none','FaceAlpha','0.1')
    end
end

%% Correction through rotation

%rotation matrix
ang=0.8036; %angle in radians, just eyeballed, and in one dimension
rm=...
    [1,0,0,0;...
    0,cos(ang),0,sin(ang);...
    0,0,1,0;...
    0,-sin(ang),0,cos(ang)]; 

%apply rotation
lsri_r=lsri(:,:)'*rm;

figure, hold on, axis equal, 
% %xlim([0 1]), ylim([-1 1]), zlim([0,2])
scatter3(lsri(1,:),lsri(2,:),lsri(4,:),[],reshape(plt_RGB,[3,220])','v','filled')
scatter3(lsri_r(:,1),lsri_r(:,2),lsri_r(:,4),[],reshape(plt_RGB,[3,220])','^','filled')

legend({'Original','Rotated'},'Location','best')
xlabel('l'),ylabel('s2'),zlabel('m2'); %l stays the same

%% Correction through shift

lsri_s = lsri; %shifted

lsri_s(4,:) = lsri(4,:)-0.27;
lsri_s(2,:) = lsri(2,:)-lsri_s(4,:);

figure, hold on, axis equal, 
scatter3(lsri(1,:),lsri(2,:),lsri(4,:),[],reshape(plt_RGB,[3,220])','v','filled')
scatter3(lsri_s(1,:),lsri_s(2,:),lsri_s(4,:),[],reshape(plt_RGB,[3,220])','^','filled')

legend({'Original','Shifted'},'Location','best')
xlabel('l'),ylabel('s2'),zlabel('m2');

