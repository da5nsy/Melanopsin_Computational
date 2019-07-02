% Prep for ICVS 2019

% Amended from melcomp_6b and melcomp_9

clear, clc, close all

set(0,'defaultAxesFontName', 'Courier')
plot_where = [100,100];
%plot_where = [1950,500];
plot_size  = [1505,727];
%plot_size  = [1100,500];

min_l_scale = 0.6;
max_l_scale = 0.9;
min_s_scale = 0;
max_s_scale = 0.1;

mktrns = 0.3; %marker transparency

base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\ICVSgifs\';
saveGifs = 0;

%% Prepare figure

f = figure('Position',[plot_where plot_size],...
    'defaultLineLineWidth',2,...
    'defaultAxesFontSize',12,...
    'WindowStyle','docked'); 
% f = figure('Position',[plot_where plot_size],...
%     'defaultLineLineWidth',2,...
%     'defaultAxesFontSize',12,...
%     'color','white'); 
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

%% Plot mean SPD

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
%xticks([0.5 1])
%yticks([0 1])
xlabel('{\itl}_{MB,10}');
ylabel('{\its}_{MB,10}');

%%
hold on
rng(1); n_ill = randi(size(T_SPD,2)); %pick a random spectrum (55th under current settings, changes if you downsample the illuminant data, for example)
scatter(lsri(1,:,n_ill),lsri(2,:,n_ill),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)

%% Zoom

interval = 50;
axlim = xlim;
bxlim = [min_l_scale, max_l_scale];
aylim = ylim;
bylim = [min_s_scale, max_s_scale];

for i=1:2
    abxlim(i,:) = linspace(axlim(i),bxlim(i),interval);
    abylim(i,:) = linspace(aylim(i),bylim(i),interval);
end

for i = 1:interval
    xlim(abxlim(:,i))
    ylim(abylim(:,i))
    drawnow
    
    if saveGifs
        frame = getframe(f);
        im{i} = frame2im(frame);
        [A,map] = rgb2ind(im{i},256);
        if i == 1
            imwrite(A,map,[base,'test','.gif'],'gif','LoopCount',0,'DelayTime',0.001);
        else
            imwrite(A,map,[base,'test','.gif'],'gif','WriteMode','append','DelayTime',0.001);
        end
    end
end

% The following doesn't work currently because I forgot that the refs were
% already a subsample of the Vrhel set (vrhel_nat_extended)

% %% add labels
% load sur_vrhel.mat
% labels_vrhel(137).label = 'peach skin -- yellow'; %correct typo
% 
% Vrhel_nat_extended_refs=[1:44,65,69,118:154];
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
    %drawnow
    % save frame of gif
end

%% Shift back to ideal

interval = 20;
a = [lsri(1,:);lsri(2,:)];
b = repmat([lsri(1,:,n_ill);lsri(2,:,n_ill)],1,size(T_SPD,2));
ab = zeros([size(a,1),size(a,2),interval]);
for i = 1:size(a,2)
    ab(1,i,:) = linspace(a(1,i),b(1,i),interval); 
    ab(2,i,:) = linspace(a(2,i),b(2,i),interval); 
end

axes(s(4))
for i=1:interval  
    cla
    scatter(ab(1,:,i),ab(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
    drawnow
end


%% Add noise

hypI = repmat([lsri(1,:,n_ill);lsri(2,:,n_ill)],1,1,130);

rng(1)
noise = normrnd(0,0.002,size(hypI));
noise(2,:,:) = noise(2,:,:)/2;

hypI = hypI+noise;
hypI2 = hypI(:,:);

interval = 20;
a = repmat([lsri(1,:,n_ill);lsri(2,:,n_ill)],1,size(T_SPD,2));
b = hypI2;
ab = zeros([size(a,1),size(a,2),interval]);
for i = 1:size(a,2)
    ab(1,i,:) = linspace(a(1,i),b(1,i),interval); 
    ab(2,i,:) = linspace(a(2,i),b(2,i),interval); 
end

axes(s(4))
for i=1:interval  
    cla
    scatter(ab(1,:,i),ab(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
    drawnow
end

%% Jump through the steps of KMeansMark

% Would be nice to show iterations rather than just result
cla
[KMM,km_idx] = KMeansMark(hypI);

colormap('lines')
scatter(hypI(1,:),hypI(2,:),[],km_idx,'filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)

text(0.8,0.08,['KMeansMark =', string(KMM)],'FontName','Courier','FontSize',15)

% Highlight correct vs incorrect

%% Now let's try some real algos...

%% But before we do, sqrt data 

%Replot original data
cla
scatter(s(4),lsri(1,:),lsri(2,:),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)

%(histograms?)
%s(5) = subplot(3,4,[2,3,4,6,7,8,10,11,12]);
%histogram(lsri(1,:,:))
%histogram(s(5),lsri(2,:,:), 'Orientation', 'horizontal')

%% Plot transformation
lsri_n = sqrt(lsri);

lsri_c = lsri_n(:,:); %would be nice to not do this transformation, for easier access later, but I'd need to be very careful to check that it didn't change the function of the calculation below

for i=1:size(lsri,1)
    lsri_c(i,:) = (lsri_c(i,:) - mean(lsri_c(i,:)))./std(lsri_c(i,:));
end

% figure,
% scatter(lsri_c(1,:),lsri_c(2,:),...
%     [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
% axis equal

lsri_n = reshape(lsri_c,size(lsri));

interval = 100;
axlim = xlim;
bxlim = [-4, 4]; %We lose a single point which goes out of bounds by setting it at +4
aylim = ylim;
%bylim = [-4, 4];
bylim = [0, 0.4];

abxlim = zeros(2,interval);
abylim = zeros(2,interval);
for i=1:2
    abxlim(i,:) = linspace(axlim(i),bxlim(i),interval);
    abylim(i,:) = linspace(aylim(i),bylim(i),interval);
end

a = [lsri(1,:);lsri(2,:)];
b = [lsri_n(1,:);lsri_n(2,:)];
ab = zeros([size(a,1),size(a,2),interval]);
for i = 1:size(a,2)
    ab(1,i,:) = linspace(a(1,i),b(1,i),interval); 
    ab(2,i,:) = linspace(a(2,i),b(2,i),interval); 
end

axes(s(4))
%axis auto
for i=1:interval  
    cla
    scatter(ab(1,:,i),ab(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
    %xlim(abxlim(:,i))
    ylim(abylim(:,i))
    drawnow
end

%scatter(s(4),lsri_n(1,:),lsri_n(2,:),'r','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)

% axis auto
% cleanTicks
% ylabel('sqrt({\its}_{MB,10})');

%%

sl2 = spectral_locus;
sl2(2,:) = sl2(2,:).^(1/2.5);
scatter(s(4),sl2(1,:),sl2(2,:),[],RGB','filled');

%% Grey world correction

% Highlight data for a single illuminant
scatter(lsri_n(1,:,n_ill),lsri_n(2,:,n_ill),'r','filled')

% Calculate mean
scatter(mean(lsri_n(1,:,n_ill)),mean(lsri_n(2,:,n_ill)),'r*')

% Calculate means for all the illuminants
lsri_m2 = mean(lsri_n,2);
scatter(squeeze(lsri_m2(1,:,:)),squeeze(lsri_m2(2,:,:)),'r*')

%%
% Drag them back together
% (change axis labels)

corrector = lsri_m2-lsri_m2(:,:,1);
lsri_mc = lsri_m2 - corrector;
lsri_c = lsri_n - corrector;
cla
scatter(lsri_c(1,:),lsri_c(2,:),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
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
    scatter3(s(4),lsri_n(1,:,i),lsri_n(2,:,i),lsri_n(4,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
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