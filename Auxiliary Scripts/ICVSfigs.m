% Prep for ICVS 2019

% Amended from melcomp_6b and melcomp_9

clear, clc, close all

set(0,'defaultAxesFontName', 'Courier')
%plot_where = [100,100];
plot_where = [1950,500];
%plot_size  = [1505,727];
plot_size  = [1100,500];

min_l_scale = 0.6;
max_l_scale = 0.9;
max_s_scale = 0.1;

mktrns = 0.3; %marker transparency

%% Prepare figure

% figure('Position',[plot_where plot_size],...
%     'defaultLineLineWidth',2,...
%     'defaultAxesFontSize',12,...
%     'WindowStyle','docked') 
figure('Position',[plot_where plot_size],...
    'defaultLineLineWidth',2,...
    'defaultAxesFontSize',12) 
hold on

s(1) = subplot(3,4,1);
s(2) = subplot(3,4,5);
s(3) = subplot(3,4,9);
s(4) = subplot(3,4,[2,3,4,6,7,8,10,11,12]);

%% Load data 

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_extended',...
    'SSF','SS10',...
    'mel_offset',0,...
    'lum','CIE_10');

%% Plot SPDs

plot(s(1),SToWls(S_sh),mean(T_SPD,2))

%% Compute colorimetry

%First compute colorimetry under a mean illuminant, to work out which
%reflectances to exclude
[~, lsri_avIll] = melcomp_colorimetry(mean(T_SPD,2), T_SRF, T_SSF, T_lum, S_sh);
% Then excludes anything closer than X to another point
exclRef = excludeSimilarRefs(lsri_avIll,0.003);
T_SRF_reduced = T_SRF(:,~exclRef);

[~, lsri] = melcomp_colorimetry(T_SPD, T_SRF_reduced, T_SSF, T_lum, S_sh);

%% Plot SRFs

plot(s(2),SToWls(S_sh),T_SRF_reduced)

%% Plot SSFs

plot(s(3),SToWls(S_sh),T_SSF(:,[1:3,5]))

%% Plot MB chromaticity diagram

% compute chromaticities of points on spectral locus
spectral_locus = LMSToMacBoyn(T_SSF(:,1:3)',T_SSF(:,1:3)',T_lum');

% compute display colours for points on spectral locus
load T_xyz1931.mat
T_xyz1931 = SplineCmf(S_xyz1931,T_xyz1931,S_sh);
RGB = XYZToSRGBPrimary(T_xyz1931);
RGB(RGB<0) = 0;
RGB(RGB>1) = 1;

scatter(s(4),spectral_locus(1,:),spectral_locus(2,:),[],RGB','filled')
%scatter(s(4),spectral_locus(1,:),spectral_locus(2,:),'k','filled') %black version
xlim([0.5 1])
ylim([0 1])
xticks([0.5 1])
yticks([0 1])
xlabel('{\itl}_{MB,10}');
ylabel('{\its}_{MB,10}');

%%
hold on
rng(1); n_ill = randi(size(T_SPD,2)); %pick a random spectrum (55th under current settings, changes if you downsample the illuminant data, for example)
scatter(lsri(1,:,n_ill),lsri(2,:,n_ill),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)

xlim([min_l_scale max_l_scale])
ylim([0 max_s_scale])
xticks([min_l_scale max_l_scale])
yticks([0 max_s_scale])

% %% add labels
% load sur_vrhel.mat
% labels_vrhel(137).label = 'peach skin -- yellow'; %correct typo
% 
% temp = 1:size(T_SRF,2);
% refs = temp(~exclRef);
% 
% for i=1:length(refs)
%     text(lsri(1,i,n_ill)+0.005,lsri(2,i,n_ill)+0.00015,labels_vrhel(refs(i)).label,'Rotation',5,'FontName','Courier')
% end

%% Add more illums

axes(s(1)), hold on

for i=1:size(T_SPD,2)
    plot(s(1),SToWls(S_sh),T_SPD(:,i))
    scatter(s(4),lsri(1,:,i),lsri(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
    drawnow
    % save frame of gif
end

%% Shift back to ideal

axes(s(4))
cla
scatter(lsri(1,:,n_ill),lsri(2,:,n_ill),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)


%% Add noise

hypI = repmat([lsri(1,:,n_ill);lsri(2,:,n_ill)],1,1,130);

rng(1)
noise = normrnd(0,0.002,size(hypI));
noise(2,:,:) = noise(2,:,:)/2;

hypI = hypI+noise;

scatter(hypI(1,:),hypI(2,:),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)

%% Jump through the steps of KMeansMark

% Would be nice to show iterations rather than just result
cla
[KMM,km_idx] = KMeansMark(hypI);

colormap('lines')
scatter(hypI(1,:),hypI(2,:),[],km_idx,'filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)

text(0.8,0.08,['KMeansMark =', string(KMM)],'FontName','Courier','FontSize',15)

% Highlight correct vs incorrect

%% Now let's try some real algos...

%% But before we do, log data 



cla

for i=1:size(T_SPD,2)
    plot(s(1),SToWls(S_sh),T_SPD(:,i))
    scatter(s(4),lsri(1,:,i),lsri(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
end

%(histograms?)
%s(5) = subplot(3,4,[2,3,4,6,7,8,10,11,12]);
%histogram(lsri(1,:,:))
%histogram(s(5),lsri(2,:,:), 'Orientation', 'horizontal')



%%
lsri = log(lsri);

lsri_c = lsri(:,:); %would be nice to not do this transformation, for easier access later, but I'd need to be very careful to check that it didn't change the function of the calculation below

for i=1:size(lsri,1)
    lsri_c(i,:) = (lsri_c(i,:) - mean(lsri_c(i,:)))./std(lsri_c(i,:));
end

% figure,
% scatter(lsri_c(1,:),lsri_c(2,:),...
%     [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
% axis equal

lsri = reshape(lsri_c,size(lsri));

cla

for i=1:size(T_SPD,2)
    plot(s(1),SToWls(S_sh),T_SPD(:,i))
    scatter(s(4),lsri(1,:,i),lsri(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
end

axis auto
cleanTicks
xlabel('norm(log({\itl}_{MB,10}))');
ylabel('norm(log({\its}_{MB,10}))');

%% Grey world correction

% Highlight data for a single illuminant
scatter(lsri(1,:,n_ill),lsri(2,:,n_ill),'r','filled')

% Calculate mean
scatter(mean(lsri(1,:,n_ill)),mean(lsri(2,:,n_ill)),'r*')

% Calculate means for all the illuminants
lsri_m2 = mean(lsri,2);
scatter(squeeze(lsri_m2(1,:,:)),squeeze(lsri_m2(2,:,:)),'r*')

%%
% Drag them back together
% (change axis labels)

cla

corrector = lsri_m2-lsri_m2(:,:,1);
lsri_mc = lsri_m2 - corrector;
lsri_c = lsri - corrector;
cla
scatter(spectral_locus(1,:),spectral_locus(2,:),'k','filled')
scatter(lsri_c(1,:),lsri_c(2,:),[],[0.5,0.5,0.5],'filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
scatter(lsri_mc(1,1,1),lsri_mc(2,1,1),'r*')
xlabel('{\it }_{ }');
ylabel('{\it }_{ }');

% Assessment

% Show k-means algo running
% Show marking of k-means

%% Same for bright-is white

%% Now for melanopsin based algo

% Handle data is exactly the same way

cla

for i=1:size(T_SPD,2)
    plot(s(1),SToWls(S_sh),T_SPD(:,i))
    scatter3(s(4),lsri(1,:,i),lsri(2,:,i),lsri(4,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
end
cleanTicks

% Shift to 3D (show sweet angle)
% Find a way to visualise the differencing

% Show the results of the 4 algos, with inset scores.

%% But what happens if we limit the number of surfaces under each illuminant

% (This is a bit like walking around with a macbeth colour checker 5 inches
% from your face)

% Then we see degraded performance from the GW and BiW algos
% This is because as we decrease the number of surfaces, we decrease the
% odds of one of these high level metrics being representative of the
% illuminant

% And so we make a prediction from this data - that the melanopsin signal
% is useful in scenes chromatic biases

%%

% But is it more valuable than another random signal?

%%

% load T_cones_ss10.mat
% load T_melanopsin.mat
% 
% %%
% 
% figure('Position',[plot_where plot_size]), hold on
% plot(SToWls(S_cones_ss10),T_cones_ss10)
% plot(SToWls(S_melanopsin),T_melanopsin)
% 
% xlim([300 900])
% yticks([min(ylim) max(ylim)])
% 
% %%
% range = -88:20:162;
% 
% for i=range
%     plot(SToWls(S_melanopsin+[i,0,0]),T_melanopsin,'k')
% end