function melcomp(PF_SPD,PF_refs,PF_obs,Z_ax)

% A fresh attempt

% TO DO
% - check that lm work out the same when I calculate them without the PTB
%   function, just for lolz (seriously - to be careful)
% - what about log signals?
% - Can the 'correction through shift' be done by division rather than subtraction?
% - Should I include Foster+'s "Levada scene"
%   (https://personalpages.manchester.ac.uk/staff/d.h.foster/Time-Lapse_HSIs/Time-Lapse_HSIs_2015.html)
%   I'd want to crop the houses out, and then we just have foliage, is this
%   neccessarily a 'salient object'?
% - Include the Foster+ 2004 images (slightly difficult because of
%   different sampling ranges/intervals)

%% Pre-flight checks
% Setting things here controls what data is used and in what way
% If we're inside a function, it is assumed that all are set, otherwise it
%   looks for them to be set manually

try
    nargin;
catch
    clear, clc, close all
    
    PF_SPD = 1;
    % 1 = CIE D series
    % 2 = Hernández-Andrés+
    
    PF_refs = 1;
    % 1 = Vhrel+ (natural only)
    % 2 = Vhrel+ (all)
    % 3 = Foster+
    
    PF_obs = 1;
    % 1 = PTB Smith-Pokorny
end

%% Load Daylight SPD
plt_SPD = 0;

if PF_SPD == 1
    D_CCT=1./linspace(1/3600,1/25000,20); %non-linear range, aiming to better reproduce observed variation
    load B_cieday
    T_SPD = GenerateCIEDay(D_CCT,[B_cieday]); %these appear to be linearly upsampled from 10nm intervals (see 'cieday investigation.m' https://github.com/da5nsy/General-Purpose-Functions/blob/3ee587429e9c4f3dd52d64acd95acf82d7e05f47/cieday%20investigation.m)
    T_SPD = (T_SPD./repmat(max(T_SPD),81,1)); %normalise
    S_SPD = S_cieday;
elseif PF_SPD == 2
    load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
    % http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
    % From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
    %       spectroradiometric characteristics of narrow-field-of-view
    %       clear skylight in Granada, Spain" (2001)
    T_SPD=final; clear final
    S_SPD=[300,5,161];
    
    if plt_SPD
        figure, hold on;        
        plot(SToWls(S_SPD),T_SPD)
        %drawnow, pause(0.3)        
        xlabel('Wavelength (nm)')
        ylabel('Spectral Power Distribution (W m ^-^2 nm^-^1)')
        xlim([S_SPD(1),max(SToWls(S_SPD))]), ylim([0 max(T_SPD(:))]);
    end
else 
    error('SPD selection failed')
end

% Load Reflectances
if or((PF_refs == 1),(PF_refs == 2))
    load sur_vrhel
    refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
    T_refs = sur_vrhel';
    if PF_refs == 1
        T_refs = sur_vrhel(:,refs)';
        T_refs = T_refs([2,1,3,6,5,4,8,9,11,10,7],:); %change order for clearer plotting visualisation
    end
    S_refs = S_vrhel;
    clear sur_vrhel refs S_vrhel 
elseif PF_refs == 3
    base = 'C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Foster Images\';
    for i=1:4 %2002 images
        ims(i)=load([base, '2002\scene',num2str(i),'.mat']); %imageS
    end
%     %2004 images
%     ims(5)=load([base,'2004\scene1\ref_crown3bb_reg1.mat']);
%     ims(6)=load([base,'2004\scene2\ref_ruivaes1bb_reg1.mat']);
%     ims(7)=load([base,'2004\scene3\ref_mosteiro4bb_reg1.mat']);
%     ims(8)=load([base,'2004\scene4\ref_cyflower1bb_reg1.mat']);
%     ims(9)=load([base,'2004\scene5\ref_cbrufefields1bb_reg1.mat']);
    
    [r, c, w] = size(ims(1).reflectances);
    T_refs = reshape(ims(1).reflectances, r*c, w);
    for i=2:4%length(ims)
        [r, c, w] = size(ims(i).reflectances);
        T_refs = [T_refs; reshape(ims(i).reflectances, r*c, w)];
    end
    S_refs=[410,10,31];
else
    error('refs selection failed')
end


% Load Observer
if PF_obs == 1
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
S_mel = S_melanopsin - [10, 0, 0]; clear S_melanopsin T_melanopsin

% For messing with hypothetical spectral sensitivity of melanopsin
%S_mel(1)=S_melanopsin(1)+30;

% Pull all observer elements together

S_sh = [max([S_SPD(1),S_refs(1),S_obs(1)]),max([S_SPD(2),S_refs(2),S_obs(2)]),min([S_SPD(3),S_refs(3),S_obs(3)])]; %S_shared: work out what the lowest common denominator for the range/interval of the data is

%reduce all data down to the common range/interval
T_SPD  = SplineSpd(S_SPD,T_SPD,S_sh,1); % extend == 1: Cubic spline, extends with last value in that direction
if or((PF_refs == 1),(PF_refs == 2))
    T_refs = SplineSrf(S_refs,T_refs',S_sh,1);
elseif PF_refs == 3 %It would theoretically be fine to run the Foster data through the above, I'm certain nothing would change, but it takes an incredibly long time (to do nothing, since the Foster data is probably always going to be the lowest resolution thing that is being handled, thus S_Sh should == S_refs)
    T_refs = T_refs';
    T_refs = T_refs(:,1:10000:end); %temporary(?) downsampling for speed in testing
    if PF_SPD == 2
        T_SPD = T_SPD(:,1:30:end);  %temporary(?) downsampling of SPD data for use with Foster+ refs
    end
end
T_obs  = SplineCmf(S_obs,T_obs,S_sh,1)';
T_rods = SplineCmf(S_rods,T_rods,S_sh)';
T_mel  = SplineCmf(S_mel,T_mel,S_sh)';
[S_SPD, S_refs, S_obs, S_rods, S_mel] = deal(S_sh);

% combine sensitivity data
T_LMSRI=[T_obs,T_rods,T_mel];
S_LMSRI=S_sh;

%% Compute colorimetry

plt_fig       = 1; % 0 = off, 1 = on
plt_locus     = 1; % plot spectral locus in the MB diagram, 0 = off, 1 = on
plt_real_cols = 0; % 0 = colormap, 1 = colours calculated from SPD and refs
plt_lines     = 0; %plot lines on the graph connecting points (good for when using low numbers of reflectances for seeing shapes through the data) 

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_refs.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_LMSRI'*T_rad(:,:,i);
    lsri(1:2,:,i) = LMSToMacBoyn(LMSRI(1:3,:,i));    
    lsri(3,:,i)   = LMSRI(4,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); 
    lsri(4,:,i)   = LMSRI(5,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); 
    % used the same scalars for luminance as are in the LMSToMAcBoyn
    %   function - which relate to the Smith-Pokorny fundamentals
end

if plt_real_cols
    % Compute colorimetry (just for display)
    load T_xyz1931.mat
    T_xyz1931=SplineCmf(S_xyz1931,T_xyz1931,S_refs)';
    for i=[11, 1:size(T_SPD,2)] %starts with 11 (6551K, arbitrary), so that a fixed white is already calculated in time for the first 'plt_RGB' line
        pltc_whiteXYZ(:,i) = T_xyz1931' * T_SPD(:,i);
        pltc_XYZ(:,:,i)    = T_xyz1931' * T_rad(:,:,i);
        pltc_Lab(:,:,i)    = XYZToLab(squeeze(pltc_XYZ(:,:,i)),pltc_whiteXYZ(:,i));
        pltc_RGB(:,:,i)    = XYZToSRGBPrimary(LabToXYZ(pltc_Lab(:,:,i),pltc_whiteXYZ(:,11))); %Using fixed, arbitrary (mid-range), white.
        %pltc_RGB(:,:,i)    = pltc_RGB(:,:,i)/max(max(pltc_RGB(:,:,i)));
    end
    pltc_RGB = pltc_RGB/max(pltc_RGB(:));
end
pltc_alt = repmat(jet(size(T_refs,2))',1,1,size(T_SPD,2)); %despite the effort gone through above to calculate the actual colours, this seems more useful for differentiating different reflectances from eachother
rng(7);
pltc_alt=pltc_alt(:,randperm(size(T_refs,2)),:); %this particular random permutation seems to generate colours in an order which means that when plotted (Hernández-Andrés+, Vrhel+) the different refs are most easily distinguishable.


try
    nargin == 4;
catch
    Z_ax = 9;
    disp('default: Z-axis is ''i''')
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
    
    if plt_real_cols == 1
        scatter3(lsri(1,:),lsri(2,:),t_Z(1,:),[],pltc_RGB(:,:)','v','filled') %with colours of objects
        hold on
    else
        scatter3(lsri(1,:),lsri(2,:),t_Z(1,:),[],pltc_alt(:,:)','v','filled') %with arbitrary colours
        hold on
    end
    if plt_lines
        plot3(lsri(1,:),lsri(2,:),t_Z(1,:),'Color',[0,0,0,0.2]) %transulent lines
    end
    zlabel(plt_lbls{Z_ax})

    if plt_locus
        MB_locus = LMSToMacBoyn(T_obs');
        %plot(MB_locus(1,:),MB_locus(2,:))
        fill([MB_locus(1,5:65),MB_locus(1,5)],[MB_locus(2,5:65),MB_locus(2,5)],'k','LineStyle','none','FaceAlpha','0.1')
    end
    
    grid on
    %axis equal
    xlim([0 1]), ylim([0 1])
    xlabel('l'),ylabel('s'), 
    %title(plt_lbls{Z_ax})
    %view(3) %view(188,46)
end

try
    nargin %ends function here
    return
catch
end

%% Correction through rotation

plt_CTR = 1;

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
    figure, hold on, axis equal, grid on
    % %xlim([0 1]), ylim([-1 1]), zlim([0,2])
    scatter3(lsri(1,:),lsri(2,:),lsri(4,:),[],pltc_alt(:,:)','v','filled')
    scatter3(lsri_r(:,1),lsri_r(:,2),lsri_r(:,4),[],pltc_alt(:,:)','^','filled')
    
    legend({'Original','Rotated'},'Location','best')
    xlabel('l'),ylabel('s2'),zlabel('i2'); %l stays the same
    
    view(90,0)
end

%% Correction through shift

plt_CTS = 1;

lsri_s = lsri; %shifted

lsri_s(2,:) = lsri(2,:)-(lsri_s(4,:)-0.27);

if plt_CTS
    figure, hold on, axis equal, grid on
    scatter3(lsri(1,:),lsri(2,:),lsri(4,:),[],pltc_alt(:,:)','v','filled')
    scatter3(lsri_s(1,:),lsri_s(2,:),lsri_s(4,:),[],pltc_alt(:,:)','^','filled')
    
    legend({'Original','Shifted'},'Location','best')
    xlabel('l'),ylabel('s2'),zlabel('i2');
    
    view(90,0)
end

end