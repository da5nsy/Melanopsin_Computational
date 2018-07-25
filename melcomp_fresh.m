function melcomp_fresh(Z_ax)

% A fresh attempt

% TO DO
% - check that lm work out the same when I calculate them without the PTB
%   function, just for lolz (seriously - to be careful)
% - what about log signals?
% - Can the 'correction through shift' be done by division rather than subtraction?


%% Pre-flight checks
% Setting things here controls what data is used and in what way

%clear, clc, close all

PF_SPD = 1;
% 1 = CIE D series
% 2 = Hern�ndez-Andr�s+

PF_refs = 1;
% 1 = Vhrel+
% 2 = Ennis+
% 3 = Foster+

PF_obs = 1;
% 1 = PTB Smith-Pokorny

%% Compute Daylight SPD
if PF_SPD == 1
    D_CCT=1./linspace(1/3600,1/25000,20); %non-linear range, aiming to better reproduce observed variation
    load B_cieday
    T_SPD = GenerateCIEDay(D_CCT,[B_cieday]); %these appear to be linearly upsampled from 10nm intervals (see 'cieday investigation.m' https://github.com/da5nsy/General-Purpose-Functions/blob/3ee587429e9c4f3dd52d64acd95acf82d7e05f47/cieday%20investigation.m)
    T_SPD = (T_SPD./repmat(max(T_SPD),81,1)); %normalise
    S_SPD = S_cieday;
elseif PF_SPD == 2
    error('Danny still needs to write the section that pulled in the Hern�ndez-Andr�s+ data')
else 
    error('SPD selection failed')
end

% Load Reflectances
if PF_refs == 1
    load sur_vrhel
    refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
    %T_vrhel = sur_vrhel';
    T_refs = sur_vrhel(:,refs)';
    T_refs = T_refs([2,1,3,6,5,4,8,9,11,10,7],:); %change order for clearer plotting visualisation
    S_refs = S_vrhel;
    clear sur_vrhel refs S_vrhel 
elseif PF_SPD == 2
    error('Danny still needs to write the section that pulled in the Ennis+ data')
elseif PF_SPD == 2
    error('Danny still needs to write the section that pulled in the Foster+ data')
else
    error('refs selection failed')
end


% Observer
if PF_SPD == 1
    % Smith-Pokorny, for use with MacLeod Boynton diagram
    load T_cones_sp
    T_obs = T_cones_sp;
    S_obs = S_cones_sp;
    clear T_cones_sp S_cones_sp
else
    error('obs selection failed')
end
load T_rods
load T_melanopsin
T_mel = SplineCmf(S_melanopsin,T_melanopsin, S_melanopsin - [10, 0, 0],1); %Increasing the range of this function in case it ends up limiting the range of S_sh, and shorten variable names
S_mel = S_melanopsin - [10, 0, 0];

% For messing with hypothetical spectral sensitivity of melanopsin
%S_mel(1)=S_melanopsin(1)+30;

%% Pull them all together

S_sh = [max([S_SPD(1),S_refs(1),S_obs(1)]),max([S_SPD(2),S_refs(2),S_obs(2)]),min([S_SPD(3),S_refs(3),S_obs(3)])]; %S_shared: work out what the lowest common denominator for the range/interval of the data is

%reduce all data down to the common range/interval
T_SPD  = SplineSpd(S_SPD,T_SPD,S_sh,1); % extend == 1: Cubic spline, extends with last value in that direction
T_refs = SplineSrf(S_refs,T_refs',S_sh,1);
T_obs  = SplineCmf(S_obs,T_obs,S_sh,1)';
T_rods = SplineCmf(S_rods,T_rods,S_sh)';
T_mel  = SplineCmf(S_mel,T_mel,S_sh)';
[S_SPD, S_refs, S_obs, S_rods, S_mel] = deal(S_sh);

% combine sensitivity data
T_LMSRI=[T_obs,T_rods,T_mel];
S_LMSRI=S_sh;

%% Combine

plt_fig     = 1;
plt_locus   = 0;
plt_ps      = 0.1; %plot pause

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_refs.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_LMSRI'*T_rad(:,:,i);
    lsri(1:2,:,i) = LMSToMacBoyn(LMSRI(1:3,:,i));    
    lsri(3,:,i)   = LMSRI(4,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); 
    lsri(4,:,i)   = LMSRI(5,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); 
    % used the same scalars for luminance as are in the LMSToMAcBoyn
    %   function - which relate to the Smith-Pokorny fundamentals
end

% Compute colorimetry (just for display)
load T_xyz1931.mat
T_xyz1931=SplineCmf(S_xyz1931,T_xyz1931,S_refs)';
for i=[11, 1:size(T_SPD,2)] %starts with 11 (6551K, arbitrary), so that a fixed white is already calculated in time for the first 'plt_RGB' line
    plt_whiteXYZ(:,i) = T_xyz1931' * T_SPD(:,i);
    plt_XYZ(:,:,i)    = T_xyz1931' * T_rad(:,:,i);
    plt_Lab(:,:,i)    = XYZToLab(squeeze(plt_XYZ(:,:,i)),plt_whiteXYZ(:,i));    
    plt_RGB(:,:,i)    = XYZToSRGBPrimary(LabToXYZ(plt_Lab(:,:,i),plt_whiteXYZ(:,11))); %Using fixed, arbitrary (mid-range), white.
    plt_RGB(:,:,i)    = plt_RGB(:,:,i)/max(max(plt_RGB(:,:,i)));
end

if ~exist('Z_ax','var') %if Z_ax isn't already defined, i.e. if we are not inside a function
    Z_ax = 9; %variable to visually test different hypothesese, 9 is default
end

plt_lbls{1}  = 'L'; %writing out this way so that there's a quick reference as to which value of Z_ax does what
plt_lbls{2}  = 'M';
plt_lbls{3}  = 'S';
plt_lbls{4}  = 'R';
plt_lbls{5}  = 'I';
plt_lbls{6}  = 'l';
plt_lbls{7}  = 's';
plt_lbls{8}  = 'r';
plt_lbls{9}  = 'i';
plt_lbls{10} = 'L+M';
plt_lbls{11} = '(0.6373*L)+(0.3924*M)';
plt_lbls{12} = 'r + i';

if plt_fig
    figure, hold on, grid on
    %axis equal
    xlim([0 1]), ylim([0 1])
    xlabel('l'),ylabel('s'), title(plt_lbls{Z_ax})
    %view(3) %view(188,46)
    
    if ismember(Z_ax,1:5)
        t_Z = LMSRI(Z_ax,:,:); %temp Z
    elseif ismember(Z_ax,6:9)
        t_Z = lsri(Z_ax-5,:,:);
    elseif Z_ax == 10
        t_Z = LMSRI(1,:,:)+LMSRI(2,:,:);
    elseif Z_ax == 11
        t_Z = 0.6373*LMSRI(1,:,:)+0.3924*LMSRI(2,:,:);
    elseif Z_ax == 12
        t_Z = lsri(3,:,:)+lsri(4,:,:);
    end
    
    for i=1:size(T_SPD,2)
        plot3(lsri(1,:,i),lsri(2,:,i),t_Z(1,:,i),'k')
        scatter3(lsri(1,:,i),lsri(2,:,i),t_Z(1,:,i),[],plt_RGB(:,:,i)','v','filled')
        zlabel(plt_lbls{Z_ax}),
        view(i*5,i)
        pause(plt_ps),drawnow
    end
    
    if plt_locus
        MB_locus=LMSToMacBoyn(T_obs');
        %plot(MB_locus(1,:),MB_locus(2,:))
        fill([MB_locus(1,5:65),MB_locus(1,5)],[MB_locus(2,5:65),MB_locus(2,5)],'k','LineStyle','none','FaceAlpha','0.1')
    end
end

%% Correction through rotation

plt_CTR = 0;

%rotation matrix
ang=0.8036; %angle in radians, just eyeballed, and in one dimension
rm=...
    [1,0,0,0;...
    0,cos(ang),0,sin(ang);...
    0,0,1,0;...
    0,-sin(ang),0,cos(ang)]; 

%apply rotation
lsri_r=lsri(:,:)'*rm;

if plt_CTR
    figure, hold on, axis equal,
    % %xlim([0 1]), ylim([-1 1]), zlim([0,2])
    scatter3(lsri(1,:),lsri(2,:),lsri(4,:),[],reshape(plt_RGB,[3,220])','v','filled')
    scatter3(lsri_r(:,1),lsri_r(:,2),lsri_r(:,4),[],reshape(plt_RGB,[3,220])','^','filled')
    
    legend({'Original','Rotated'},'Location','best')
    xlabel('l'),ylabel('s2'),zlabel('i2'); %l stays the same
end

%% Correction through shift

plt_CTS = 0;

lsri_s = lsri; %shifted

lsri_s(4,:) = lsri(4,:)-0.27;
lsri_s(2,:) = lsri(2,:)-lsri_s(4,:);

if plt_CTS
    figure, hold on, axis equal,
    scatter3(lsri(1,:),lsri(2,:),lsri(4,:),[],reshape(plt_RGB,[3,220])','v','filled')
    scatter3(lsri_s(1,:),lsri_s(2,:),lsri_s(4,:),[],reshape(plt_RGB,[3,220])','^','filled')
    
    legend({'Original','Shifted'},'Location','best')
    xlabel('l'),ylabel('s2'),zlabel('i2');
end

end