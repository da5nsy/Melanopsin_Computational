% Prep for ICVS 2019

% Amended from melcomp_6b and melcomp_9

clear, clc, close all

set(0,'defaultAxesFontName', 'Courier')

min_l_scale = 0.6;
max_l_scale = 0.9;
min_s_scale = 0;
max_s_scale = 0.1;

mktrns = 0.3; %marker transparency

base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\ICVSfigs\';
saveFigs = 1;
if saveFigs
    warning('Save figs is on - you sure? This will overwrite.')
end

demo_ill = 15;

%% Prepare figure
%[1950,500]

f = figure('Position',[[100,100], [1505,727]],...
    'defaultLineLineWidth',2,...
    'defaultAxesFontSize',12,...
    'color','white'); 
hold on
clf

%% Load data 

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_extended',...
    'SSF','SS10',...
    'mel_offset',0,...
    'lum','CIE_10');

%% Plot SPD

s(1) = subplot(3,4,1);
plot(s(1),SToWls(S_sh),T_SPD(:,demo_ill)) %random number just to get one that looks good on the axes
xlim([390 730])
ylim([0 1.7])
yticks([0 1])
xlabel('Wavelength (nm)')
ylabel('Spec. Irradiance')
if saveFigs
    print([base,'1a_SPD.png'],'-dpng','-r0')
end

%% Compute colorimetry

%First compute colorimetry under a mean illuminant, to work out which
%reflectances to exclude
[~, lsri_avIll] = melcomp_colorimetry(mean(T_SPD,2), T_SRF, T_SSF, T_lum, S_sh);
% Then excludes anything closer than X to another point
exclRef = excludeSimilarRefs(lsri_avIll,0.003);
T_SRF_reduced = T_SRF(:,~exclRef);

[~, lsri] = melcomp_colorimetry(T_SPD, T_SRF_reduced, T_SSF, T_lum, S_sh);

%% Plot SRFs
s(2) = subplot(3,4,5);
plot(SToWls(S_sh),T_SRF_reduced)
xlim([390 730])
xlabel('Wavelength (nm)')
ylabel('Reflectance')
yticks(ylim)
if saveFigs
    print([base,'1b_SRFs.png'],'-dpng','-r0')
end

%% Plot SSFs
s(3) = subplot(3,4,9);
plot(SToWls(S_sh),T_SSF(:,1:3))
xlim([390 730])
xlabel('Wavelength (nm)')
ylabel('Norm. Spec. Sens.')
yticks(ylim)
if saveFigs
    print([base,'1c_SSFs.png'],'-dpng','-r0')
end

%% Plot MB chromaticity diagram

% compute chromaticities of points on spectral locus
spectral_locus = LMSToMacBoyn(T_SSF(:,1:3)',T_SSF(:,1:3)',T_lum');

% compute display colours for points on spectral locus
load T_xyz1931.mat
T_xyz1931 = SplineCmf(S_xyz1931,T_xyz1931,S_sh);
RGB = XYZToSRGBPrimary(T_xyz1931);
RGB(RGB<0) = 0;
RGB(RGB>1) = 1;

s(4) = subplot(3,4,[2,3,4,6,7,8,10,11,12]);
scatter(s(4),spectral_locus(1,:),spectral_locus(2,:),[],RGB','filled')
%scatter(s(4),spectral_locus(1,:),spectral_locus(2,:),'k','filled') %black version
xlim([0.5 1])
ylim([0 1])
%xticks([0.5 1])
%yticks([0 1])
xlabel('{\itl}_{MB,10}');
ylabel('{\its}_{MB,10}');
if saveFigs
    print([base,'1e_MB.png'],'-dpng','-r0')
end

%%
hold on
scatter(lsri(1,:,demo_ill),lsri(2,:,demo_ill),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
if saveFigs
    print([base,'1f_MBpoints.png'],'-dpng','-r0')
end

%% Zoom

interval = 20;
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
        
    if saveFigs
        createGIF(f,base,'1g_MBzoom',i) %Outside repo: available at: https://github.com/da5nsy/General-Purpose-Functions/blob/1604ee42c2377c72c99bc30c47bc839546aa7afb/createGIF.m
    end
end

%% add labels

% % The following doesn't work currently because I forgot that the refs were
% % already a subsample of the Vrhel set (vrhel_nat_extended)
 
% load sur_vrhel.mat
% labels_vrhel(137).label = 'peach skin -- yellow'; %correct typo
% 
% Vrhel_nat_extended_refs=[1:44,65,69,118:154];
% 
% temp = 1:size(T_SRF,2);
% refs = temp(~exclRef);
% 
% for i=1:length(refs)
%     text(lsri(1,i,15)+0.005,lsri(2,i,15)+0.00015,labels_vrhel(refs(i)).label,'Rotation',5,'FontName','Courier')
% end

%% Add more illums

axes(s(1)), hold on

for i=1:size(T_SPD,2)
    plot(s(1),SToWls(S_sh),T_SPD(:,i))
    scatter(s(4),lsri(1,:,i),lsri(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
    if saveFigs
        createGIF(f,base,'1h_MBaddIllums',i)
    end
end

%% Shift back to ideal

interval = 20;
a = [lsri(1,:);lsri(2,:)];
b = repmat([lsri(1,:,demo_ill);lsri(2,:,demo_ill)],1,size(T_SPD,2));
ab = zeros([size(a,1),size(a,2),interval]);
for i = 1:size(a,2)
    ab(1,i,:) = linspace(a(1,i),b(1,i),interval); 
    ab(2,i,:) = linspace(a(2,i),b(2,i),interval); 
end

axes(s(4))
for i=1:interval  
    cla
    scatter(ab(1,:,i),ab(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
    if saveFigs
        createGIF(f,base,'1i_MBshiftBackToIdeal',i)
    end
end


%% Add noise

hypI = repmat([lsri(1,:,demo_ill);lsri(2,:,demo_ill)],1,1,130);

rng(1)
noise = normrnd(0,0.002,size(hypI));
noise(2,:,:) = noise(2,:,:)/2;

hypI = hypI+noise;
hypI2 = hypI(:,:);

interval = 20;
a = repmat([lsri(1,:,demo_ill);lsri(2,:,demo_ill)],1,size(T_SPD,2));
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
    if saveFigs
        createGIF(f,base,'1j_MBaddNoise',i)
    end
end

if saveFigs
    print([base,'1j2_MBNoise.png'],'-dpng','-r0')
end

%% Jump through the steps of KMeansMark

% Would be nice to show iterations rather than just result
cla
[KMM,km_idx] = KMeansMark(hypI);

colormap('lines')
scatter(hypI(1,:),hypI(2,:),[],km_idx,'filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)

text(0.8,0.08,['KMeansMark =', string(KMM)],'FontName','Courier','FontSize',15)

if saveFigs
    print([base,'1k_KMM.png'],'-dpng','-r0')
end

% Highlight correct vs incorrect

%% Now let's try some real algos...

%% But before we do, log data 

%Replot original data
cla
scatter(s(4),lsri(1,:),lsri(2,:),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)

%(histograms?)
%s(5) = subplot(3,4,[2,3,4,6,7,8,10,11,12]);
%histogram(lsri(1,:,:))
%histogram(s(5),lsri(2,:,:), 'Orientation', 'horizontal')

%% Plot transformation
lsri_n = log(lsri);

lsri_c = lsri_n(:,:); 
for i=1:size(lsri,1)
    lsri_c(i,:) = (lsri_c(i,:) - mean(lsri_c(i,:)))./std(lsri_c(i,:));
end
lsri_n = reshape(lsri_c,size(lsri));

interval = 20;
axlim = xlim;
bxlim = [-3, 5]; 
aylim = ylim;
bylim = [-3, 5];

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
for i=1:interval  
    cla
    scatter(ab(1,:,i),ab(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
    xlim(abxlim(:,i))
    ylim(abylim(:,i))
    drawnow
    if saveFigs
        createGIF(f,base,'1l_PreAlgoTransf',i)
    end
end

axis equal %!!!!!!!!!!!!!!!!!! check this later, and pull the assigned values into bxlim and bylim
xlabel('norm(log({\itl}_{MB,10}))');
ylabel('norm(log({\its}_{MB,10}))');

%% Check effect on spectral locus

% Currently not working as spectral locus isn't normalised with the signals

% sl2 = spectral_locus;
% sl2(2,:) = sqrt(sl2(2,:));
% scatter(s(4),sl2(1,:),sl2(2,:),[],RGB','filled');
% axis auto

%% Perform CC

% Lum = zeros(size(T_SRF_reduced,2),size(T_SPD,2));
% for i=1:size(lsri,3)
%     Lum(:,i) = T_lum'*(T_SRF_reduced.*T_SPD(:,i));
% end
% 
% output = performCC(lsri_n,Lum);

% output_norm(1,1:nSurf(pcSurf),:,:,pcSurf) =  output(1,1:nSurf(pcSurf),:,:,pcSurf)./std(output(1,1:nSurf(pcSurf),:,:,pcSurf));
% output_norm(2,1:nSurf(pcSurf),:,:,pcSurf) =  output(2,1:nSurf(pcSurf),:,:,pcSurf)./std(output(1,1:nSurf(pcSurf),:,:,pcSurf));

% figure,
% scatter(reshape(output(1,:,:,4),[],1),reshape(output(2,:,:,4),[],1))

%% Algo visualisation

f2 = figure('Position',[[100,100], [1505,727]],...
    'defaultLineLineWidth',2,...
    'defaultAxesFontSize',12,...
    'color','white'); 
hold on
clf

% DN
s2(1) = subplot(2,2,1);
hold on
axis equal
scatter(lsri_n(1,:),lsri_n(2,:),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)

if saveFigs
    print([base,'2a_DN.png'],'-dpng','-r0')
end

%% GW
s2(2) = subplot(2,2,2);
hold on
axis equal

scatter(lsri_n(1,:),lsri_n(2,:),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
if saveFigs
    print([base,'2b_GW0.png'],'-dpng','-r0')
end

% Highlight data for a single illuminant
scatter(s2(2),lsri_n(1,:,demo_ill),lsri_n(2,:,demo_ill),'r','filled')
if saveFigs
    print([base,'2b_GW1.png'],'-dpng','-r0')
end

% Calculate mean
scatter(s2(2),mean(lsri_n(1,:,demo_ill)),mean(lsri_n(2,:,demo_ill)),100,'b^','filled')
if saveFigs
    print([base,'2b_GW2.png'],'-dpng','-r0')
end

% Calculate means for all the illuminants
%cla
%scatter(lsri_n(1,:),lsri_n(2,:),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns) %replot original data
lsri_m2 = mean(lsri_n,2);
scatter(s2(2),squeeze(lsri_m2(1,:,:)),squeeze(lsri_m2(2,:,:)),100,'b^','filled')
if saveFigs
    print([base,'2b_GW3.png'],'-dpng','-r0')
end

lsri_c = lsri_n - lsri_m2;

interval = 20;
a = [lsri_n(1,:);lsri_n(2,:)];
b = [lsri_c(1,:);lsri_c(2,:)];
ab = zeros([size(a,1),size(a,2),interval]);
for i = 1:size(a,2)
    ab(1,i,:) = linspace(a(1,i),b(1,i),interval); 
    ab(2,i,:) = linspace(a(2,i),b(2,i),interval); 
end

c = [squeeze(lsri_m2(1,:,:))';squeeze(lsri_m2(2,:,:))'];
%d = repmat(lsri_m2(1:2,:,1),1,130);
d = zeros(size(c));
cd = zeros([size(c,1),size(c,2),interval]);
for i = 1:size(c,2)
    cd(1,i,:) = linspace(c(1,i),d(1,i),interval); 
    cd(2,i,:) = linspace(c(2,i),d(2,i),interval); 
end

cla
for i=1:interval  
    cla
    scatter(ab(1,:,i),ab(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
    scatter(cd(1,:,i),cd(2,:,i),100,'b^','filled')
    drawnow
    if saveFigs
        createGIF(f2,base,'2b_GW4',i)
    end
end

%plot([0,0],ylim,'k:')
%plot(xlim,[0,0],'k:')
if saveFigs
    print([base,'2b_GW5.png'],'-dpng','-r0')
end

%% Same for bright-is white

s2(3) = subplot(2,2,3);
hold on
axis equal

scatter(lsri_n(1,:),lsri_n(2,:),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
if saveFigs
    print([base,'2c_BiW1.png'],'-dpng','-r0')
end

% Calculate Luminance
Lum = zeros(size(T_SRF_reduced,2),size(T_SPD,2));
for i=1:size(lsri,3)
    Lum(:,i) = T_lum'*(T_SRF_reduced.*T_SPD(:,i));
end

% Highlight data for a single illuminant
scatter(s2(3),lsri_n(1,:,demo_ill),lsri_n(2,:,demo_ill),'r','filled')
if saveFigs
    print([base,'2c_BiW2.png'],'-dpng','-r0')
end

[~,maxLoc] = max(Lum);
scatter(s2(3),lsri_n(1,maxLoc(demo_ill),demo_ill),lsri_n(2,maxLoc(demo_ill),demo_ill),100,'b^','filled')
if saveFigs
    print([base,'2c_BiW3.png'],'-dpng','-r0')
end

% Calculate maxes for all the illuminants
cla
for i=1:size(T_SPD,2)
    lsri_biwm(:,1,i) = lsri_n(:,maxLoc(i),i); %lsri Bright is White max
end
scatter(lsri_n(1,:),lsri_n(2,:),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns) %replot original data
scatter(s2(3),squeeze(lsri_biwm(1,:,:)),squeeze(lsri_biwm(2,:,:)),100,'b^','filled')
if saveFigs
    print([base,'2c_BiW4.png'],'-dpng','-r0')
end

lsri_c_biw = lsri_n - lsri_biwm;

interval = 20;
a = [lsri_n(1,:);lsri_n(2,:)];
b = [lsri_c_biw(1,:);lsri_c_biw(2,:)];
ab = zeros([size(a,1),size(a,2),interval]);
for i = 1:size(a,2)
    ab(1,i,:) = linspace(a(1,i),b(1,i),interval); 
    ab(2,i,:) = linspace(a(2,i),b(2,i),interval); 
end

c = [squeeze(lsri_biwm(1,:,:))';squeeze(lsri_biwm(2,:,:))'];
%d = repmat(lsri_m2(1:2,:,1),1,130);
d = zeros(size(c));
cd = zeros([size(c,1),size(c,2),interval]);
for i = 1:size(c,2)
    cd(1,i,:) = linspace(c(1,i),d(1,i),interval); 
    cd(2,i,:) = linspace(c(2,i),d(2,i),interval); 
end

cla
for i=1:interval  
    cla
    scatter(ab(1,:,i),ab(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
    scatter(cd(1,:,i),cd(2,:,i),100,'b^','filled')
    drawnow
    if saveFigs
        createGIF(f2,base,'2c_BiW5',i)
    end
end

%plot([0,0],ylim,'k:')
%plot(xlim,[0,0],'k:')
if saveFigs
    print([base,'2b_BiW6.png'],'-dpng','-r0')
end

%% Now for melanopsin based algo

s2(4) = subplot(2,2,4);
hold on
axis equal

scatter(lsri_n(1,:),lsri_n(2,:),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
if saveFigs
    print([base,'2d_MC1.png'],'-dpng','-r0')
end

% Highlight data for a single illuminant
scatter(s2(4),lsri_n(1,:,demo_ill),lsri_n(2,:,demo_ill),'r','filled')
if saveFigs
    print([base,'2d_MC2.png'],'-dpng','-r0')
end

% Highlight data for a single illuminant
scatter(s2(4),lsri_n(1,:,demo_ill),lsri_n(2,:,demo_ill),(lsri_n(4,:,demo_ill)-min(lsri_n(4,:,demo_ill))+0.0001)*20,'b^','filled')
if saveFigs
    print([base,'2d_MC3.png'],'-dpng','-r0')
end

[sf_l,sf_s] = melcomp_6_calcsf(lsri_n, -10:0.001:10,-10:0.001:10); %calculates scaling factors

% sf_l = 0.6890;
% sf_s = -0.9300;
MC = lsri_n;
MC(1,:) = MC(1,:)+sf_l*MC(4,:);
MC(2,:) = MC(2,:)+sf_s*MC(4,:);

% cla
% scatter(s2(4),MC(1,:),MC(2,:),'k','filled');

interval = 20;
a = [lsri_n(1,:);lsri_n(2,:)];
b = MC(1:2,:);
ab = zeros([size(a,1),size(a,2),interval]);
for i = 1:size(a,2)
    ab(1,i,:) = linspace(a(1,i),b(1,i),interval); 
    ab(2,i,:) = linspace(a(2,i),b(2,i),interval); 
end

for i=1:interval  
    cla
    scatter(ab(1,:,i),ab(2,:,i),'k','filled','MarkerFaceAlpha',mktrns,'MarkerEdgeAlpha',mktrns)
    drawnow
    if saveFigs
        createGIF(f2,base,'2d_MC4',i)
    end
end

%%

jlist = 12:100;
for pcSurf = 1:length(jlist)
    nSurf(pcSurf) = round((jlist(pcSurf)/100) * size(T_SRF_reduced,2));
end

try
    load jlist.mat
catch
    [mark,output,sel_store] = melcomp_9('plt','none','pcSurfRange',jlist);
    save('jlist.mat','mark','output','sel_store')
end

counter = 1;
for j=length(jlist):-1:1
    for i=1:4
        axes(s2(i))
        cla
        %axis manual
        gscatter(reshape(output(1,1:nSurf(j),:,i,j),[],1),...
                 reshape(output(2,1:nSurf(j),:,i,j),[],1),...
            reshape(sel_store(j,1:nSurf(j),:),[],1),...
            [],[],5)
        legend off
    end
    drawnow
    if saveFigs
        if counter == 1            
            print([base,'3a_degredation.png'],'-dpng','-r0')
        end
        createGIF(f2,base,'3b_degredation',counter)
        counter = counter + 1;
    end
end


%% Summary plot

f3 = figure('Position',[[100,100], [1505,727]],...
    'defaultLineLineWidth',2,...
    'defaultAxesFontSize',12,...
    'color','white'); 
hold on

titles = {'Do nothing','Grey World','Bright-is-White','Melanopsin'};
markers = {'s:','d:','^:','o:'};
nMethods = 4;
mfc = hsv(nMethods); %marker face colour
lw = 3; %linewidth

for i=1:nMethods
    plot(jlist,mark(i,:),markers{i},'Color',mfc(i,:),'MarkerFaceColor',mfc(i,:),'linewidth',lw)
end
legend(titles,'Location','Northwest')
%plot(p.Results.pcSurfRange,1./round((p.Results.pcSurfRange/100) * size(T_SRF_reduced,2)),'k','DisplayName','Chance','linewidth',lw-1)
xlim([0 100])
ylim([0 1])

xlabel(sprintf('Percentage of surfaces used in each run (/%d)',size(T_SRF_reduced,2)))
ylabel('K-means-mark')
grid on

if saveFigs
    print([base,'4a_summary.png'],'-dpng','-r0')
end

%%

% But is it more valuable than another random signal?

f3 = figure('Position',[[100,100], [1505,727]],...
    'defaultLineLineWidth',2,...
    'defaultAxesFontSize',12,...
    'color','white'); 

melcomp_9_caller
xlim([min(xlim),670])

if saveFigs
    print([base,'5a_diffMel.png'],'-dpng','-r0')
end

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

%%

