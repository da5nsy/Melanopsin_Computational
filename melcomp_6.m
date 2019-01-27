%% Creating visual narrative for melcomp presentation

clear, clc, close all

% Load observer data

load T_cones_ss10.mat; 
T_SSF = T_cones_ss10; 
S_SSF = S_cones_ss10;
clear T_cones_ss10 S_cones_ss10

load T_rods T_rods S_rods
load T_melanopsin T_melanopsin S_melanopsin
T_mel = SplineCmf(S_melanopsin, T_melanopsin, S_melanopsin - [10, 0, 0],1); %Increasing the range of this function in case it ends up limiting the range of S_sh, and shorten variable names
S_mel = S_melanopsin - [10, 0, 0];
clear S_melanopsin T_melanopsin

% Load reflectance data

load sur_vrhel
refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
T_SRF = sur_vrhel(:,refs);
S_SRF = S_vrhel;
clear sur_vrhel S_vrhel

% Load daylight data

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)
T_SPD=final; clear final
S_SPD=[300,5,161];

%reduce all data down to the common range/interval
S_sh = [max([S_SPD(1),S_SRF(1),S_SSF(1)]),max([S_SPD(2),S_SRF(2),S_SSF(2)]),min([S_SPD(3),S_SRF(3),S_SSF(3)])]; %S_shared: work out what the lowest common denominator for the range/interval of the data is
%S_sh = [400,5,61]; % brute force way to do something like the variableweights thing

T_SPD = SplineSpd(S_SPD,T_SPD,S_sh,1); % extend == 1: Cubic spline, extends with last value in that direction
T_SRF = SplineSrf(S_SRF,T_SRF,S_sh,1);
T_SSF  = SplineCmf(S_SSF,T_SSF,S_sh,1)';
T_rods = SplineCmf(S_rods,T_rods,S_sh)'; %extend? !!!!!!!!!!
T_mel  = SplineCmf(S_mel,T_mel,S_sh)'; %extend? !!!!!!!!!!
[S_SPD, S_SRF, S_SSF, S_rods, S_mel] = deal(S_sh);

% combine sensitivity vectors
T_LMSRI=[T_SSF,T_rods,T_mel];
S_LMSRI=S_sh;

%% Plot MB chromaticity diagram

sf_10 = [0.69283932, 0.34967567, 0.05547858]; %energy 10deg from CIE 170-2:2015
T_lum = sf_10(1)*T_SSF(:,1)+sf_10(2)*T_SSF(:,2);
spectral_locus = LMSToMacBoynDG(T_SSF',T_SSF',T_lum');
load T_xyz1931.mat
T_xyz1931 = SplineCmf(S_xyz1931,T_xyz1931,S_sh);
RGB = XYZToSRGBPrimary(T_xyz1931);
RGB(RGB<0) = 0;
RGB(RGB>1) = 1;

figure('Position',[500 200 800 400]), hold on
scatter(spectral_locus(1,:),spectral_locus(2,:),[],RGB','filled')
xlim([0.5 1])
ylim([0 1])
xticks([0.5 1])
yticks([0 1])
xlabel('{\itl}_{MB}');
ylabel('{\its}_{MB}');

%save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MB.pdf")

scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
%save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBblack.pdf")

%compute colours for display
pltc_alt = repmat(jet(size(T_SRF,2))',1,1,size(T_SPD,2)); %despite the effort gone through above to calculate the actual colours, this seems more useful for differentiating different reflectances from eachother
rng(7);
pltc_alt=pltc_alt(:,randperm(size(T_SRF,2)),:); %this particular random permutation seems to generate colours in an order which means that when plotted (Hernández-Andrés+, Vrhel+) the different refs are most easily distinguishable.


%%

T_rad = zeros([S_sh(3),size(T_SRF,2),size(T_SPD,2)]);
LMSRI = zeros([size(T_LMSRI,2),size(T_SRF,2),size(T_SPD,2)]);

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_SRF.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_LMSRI'*T_rad(:,:,i);
    lsri(1:2,:,i) = LMSToMacBoynDG(LMSRI(1:3,:,i));    
%     lsri(3,:,i)   = LMSRI(4,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); 
%     lsri(4,:,i)   = LMSRI(5,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); 
end

scatter(lsri(1,:,1000),lsri(2,:,1000),[],pltc_alt(:,:,1000)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6) %1000 is just a random illuminant, no logic behind choice
%save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBsingleset_nn.pdf")

xlim([0.6 0.8])
ylim([0 0.02])
xticks([0.6 0.8])
yticks([0 0.02])

%[],pltc_alt(:,:)','v','filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4

%save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBsingleset_n.pdf")

%%
scatter(lsri(1,:),lsri(2,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
%save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBall.pdf")

cla
scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
scatter(lsri(1,:),lsri(2,:),'k','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBallgrey.pdf")


