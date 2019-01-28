%% Creating visual narrative for melcomp presentation

%Pre-flight
clear, clc, close all

plot_where = [500,200];
plot_size  = [800,400];

min_l_scale = 0.62;
max_l_scale = 0.82;
max_s_scale = 0.04;

set(0,'defaultAxesFontName', 'Courier')

%% Data
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
refs=[87, 93, 134, 137, 138, 65, 19, 24, 140, 26];
T_SRF = sur_vrhel(:,refs);
S_SRF = S_vrhel;
clear sur_vrhel S_vrhel

% Load daylight data

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hern�ndez-Andr�s, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)
T_SPD=final; clear final
T_SPD = T_SPD(:,1:20:end);
S_SPD=[300,5,161];

%% reduce all data down to the common range/interval
S_sh = [max([S_SPD(1),S_SRF(1),S_SSF(1),S_rods(1),S_mel(1)]),...
    max([S_SPD(2),S_SRF(2),S_SSF(2),S_rods(2),S_mel(2)]),...
    min([S_SPD(3),S_SRF(3),S_SSF(3),S_rods(3),S_mel(3)])]; 
%S_shared: work out what the lowest common denominator for the range/interval of the data is

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

figure('Position',[plot_where plot_size]), hold on
scatter(spectral_locus(1,:),spectral_locus(2,:),[],RGB','filled')
xlim([0.5 1])
ylim([0 1])
xticks([0.5 1])
yticks([0 1])
xlabel('{\itl}_{MB}');
ylabel('{\its}_{MB}');

save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MB.pdf")

scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBblack.pdf")

%% Compute colorimetry

T_rad = zeros([S_sh(3),size(T_SRF,2),size(T_SPD,2)]);
LMSRI = zeros([size(T_LMSRI,2),size(T_SRF,2),size(T_SPD,2)]);
lsri  = zeros([4,size(T_SRF,2),size(T_SPD,2)]); %hard coded 4
t_r   = zeros([2,size(T_SRF,2),size(T_SPD,2)]);
t_i   = zeros([2,size(T_SRF,2),size(T_SPD,2)]);

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_SRF.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_LMSRI'*T_rad(:,:,i);
    lsri(1:2,:,i) = LMSToMacBoynDG(LMSRI(1:3,:,i),T_SSF',T_lum');
    t_r(:,:,i)    = LMSToMacBoynDG(LMSRI([1,2,4],:,i),[T_SSF(:,1:2)';T_rods'],T_lum');
    t_i(:,:,i)    = LMSToMacBoynDG(LMSRI([1,2,5],:,i),[T_SSF(:,1:2)';T_mel'],T_lum');
end
lsri(3,:,:) = t_r(2,:,:); clear t_r
lsri(4,:,:) = t_i(2,:,:); clear t_i

%compute colours for display
pltc_alt = repmat(jet(size(T_SRF,2))',1,1,size(T_SPD,2)); 
rng(7); pltc_alt=pltc_alt(:,randperm(size(T_SRF,2)),:); %best combo for differentiating close chromaticities

%plot MB with points, not normalised
rng(1); n_ill = randi(size(T_SPD,2)); %pick a random spectrum (55th)
scatter(lsri(1,:,n_ill),lsri(2,:,n_ill),[],pltc_alt(:,:,n_ill)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBsingleset_nn.pdf")

%rescale diagram
xlim([min_l_scale max_l_scale])
ylim([0 max_s_scale])
xticks([min_l_scale max_l_scale])
yticks([0 max_s_scale])
save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBsingleset_n.pdf")

%add labels
labels_vrhel(137).label = 'peach skin -- yellow';
for i=1:length(refs)
    text(lsri(1,i,n_ill)+0.005,lsri(2,i,n_ill)+0.00015,labels_vrhel(refs(i)).label,'Rotation',10,'FontName','Courier')
end
save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBsingleset_n_text.pdf")

scatter(0.699237,0.025841,'rs')

%%
% all points, in colour
cla
scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
scatter(lsri(1,:),lsri(2,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBall.pdf")

%all points in grey
cla
scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
scatter(lsri(1,:),lsri(2,:),[],[0.5,0.5,0.5],'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBallgrey.pdf")

%means under each illuminant
lsri_m2 = mean(lsri,2); 
scatter(squeeze(lsri_m2(1,:,:)),squeeze(lsri_m2(2,:,:)),'r*')
save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBmeans.pdf")

%% Grey world adaptation
corrector = lsri_m2-lsri_m2(:,:,1);
lsri_mc = lsri_m2 - corrector;
lsri_c = lsri - corrector;
cla
scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
scatter(lsri_c(1,:),lsri_c(2,:),[],[0.5,0.5,0.5],'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
scatter(lsri_mc(1,1,1),lsri_mc(2,1,1),'r*')
xlabel('{\it }_{ }');
ylabel('{\it }_{ }');
save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MBc.pdf")

%% Splits

%would be nice to bring plot on piece by piece
%and to get all of this to loop instead of being seperate

figure(2,'Position',[plot_where plot_size]), hold on

s(1) = subplot(1,2,1);
hold on
s(2) = subplot(1,2,2);
hold on

xlim(s(1),[min_l_scale max_l_scale])
xticks(s(1),[min_l_scale max_l_scale])
xlabel(s(1),'{\itl}_{MB}');
yticks(s(1),[])

xlim(s(2),[0 max_s_scale])
xticks(s(2),[0 max_s_scale])
xlabel(s(2),'{\its}_{MB}');
yticks(s(2),[])

save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\split_empty.pdf")

ylabel(s(1),'{\itI}');
scatter(s(1),lsri(1,:),LMSRI(5,:),[],[0.5,0.5,0.5],'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
scatter(s(2),lsri(2,:),LMSRI(5,:),[],[0.5,0.5,0.5],'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
yticks(s(1),[min(ylim),max(ylim)])
yticks(s(2),[])

save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\split_I_grey.pdf")

%LMSRI
plotOrderNums = [5,1,2,3];
plotOrderNames = {'I','L','M','S'};
for i=1:length(plotOrderNums)
    cla(s(1))
    cla(s(2))
    
    ylabel(s(1),['{\it',plotOrderNames{i},'}']);    
    
    scatter(s(1),lsri(1,:),LMSRI(plotOrderNums(i),:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    scatter(s(2),lsri(2,:),LMSRI(plotOrderNums(i),:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    
    yticks(s(1),[min(ylim),max(ylim)])
    yticks(s(2),[])
    
    save2pdf(['C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\split_',plotOrderNames{i},'.pdf'])
end

%lsri
plotOrderNums = [1,2,4];
plotOrderNames = {'{\itl}_{MB}','{\itl}_{MB}','{\iti}_{MB}'};
plotOrderSaveNames = {'lmb','smb','imb'};
for i=1:length(plotOrderNums)
    cla(s(1))
    cla(s(2))
    
    ylabel(s(1),plotOrderNames{i});    
    
    scatter(s(1),lsri(1,:),lsri(plotOrderNums(i),:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    scatter(s(2),lsri(2,:),lsri(plotOrderNums(i),:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    
    yticks(s(1),[min(ylim),max(ylim)])
    yticks(s(2),[])
    
    save2pdf(['C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\split_',plotOrderSaveNames{i},'.pdf'])
end

%%
figure('Position',[plot_where plot_size]), hold on
xlim([min_l_scale max_l_scale])
ylim([0 max_s_scale])
xticks([min_l_scale max_l_scale])
yticks([0 max_s_scale])
xlabel('{\itl}_{MB}');
ylabel('{\its}_{MB}');

scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')

scatter(lsri(1,:),lsri(2,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)

cla
%scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')

scatter(lsri(1,:)-lsri(1,:),lsri(2,:)-lsri(1,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
axis('auto')

xticks([min_l_scale max_l_scale])
yticks([0 max_s_scale])

%save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\MB.pdf")


%%


%save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Oxford Presentation\figs\opt.pdf")






