clear, clc, close all

% A fresh attempt

% TO DO
% - check that lm work out the same when I calculate them without the PTB
%   function, just for lolz (seriously - to be careful)

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

% For messing with hypothetical spectral sensitivity of melanopsin
%S_melanopsin(1)=S_melanopsin(1)+30;

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

%% Rotation of data

lsm=[ls(1:2,:);m(:)']';

%rotation matrix



a=0.8036; %angle in radians, just eyeballed
rm=...
    [1,0,0;...
    0,cos(a),-sin(a);...
    0,sin(a),cos(a)]; 

lsm_r=lsm*rm';

figure, hold on, axis equal, 
xlim([0 1]), ylim([-1 1]), zlim([0,1])
plot3(lsm(:,1),lsm(:,2),lsm(:,3),'bo')
plot3(lsm_r(:,1),lsm_r(:,2),lsm_r(:,3),'ro')

%view(188,46)

%% %-% %%
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

ls4 = ls./[ls(1,:,:);ls(2,:,:)];
figure, hold on, axis equal, %xlim([0 1]), ylim([0 1])
for i=1:4:size(T_Dspd,1)
    plot(ls4(1,:,i),ls4(2,:,i),'o-')
    title('Yes, it really is meant to just be a single point')
end

% ls5 = ls./ls; %same as above, just being careful
% figure, hold on, axis equal, %xlim([0 1]), ylim([0 1])
% for i=1:4:size(T_Dspd,1)
%     plot(ls5(1,:,i),ls5(2,:,i),'o-')
%     title('Yes, it really is meant to just be a single point')
% end


%% calibrating by m

ls6 = ls./[m(1,:,:);m(1,:,:)];
%ls6(2,:)=(1./(ls6(2,:)+0.1245)).^(1/0.44);
figure, hold on, axis equal, %xlim([0 1]), ylim([0 1])
for i=1:4:size(T_Dspd,1)
    plot(ls6(1,:,i),ls6(2,:,i),'o-')
end
drawnow

% What would this ideally look like?

%% Trying to work out the equation

a=squeeze(ls6(:,4,:));

figure, hold on
clf, hold on
scatter(a(1,:),a(2,:))
x=1:0.1:3.1;
y=(1./(x.^0.44))-0.1245;
plot(x,y)

% y2=(1./(y+0.1245)).^(1/0.44);
% plot(x,y2)
% %%
% a2=a;
% a2(2,:)=(1./(a2(2,:)+0.1245)).^(1/0.44);
% 
% figure,
% scatter(a2(1,:),a2(2,:))


