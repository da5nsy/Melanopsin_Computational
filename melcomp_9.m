% The goal here is to compare a melanopsin based correction with other
% basic corrections (do nothing, grey world, bright-is-white/max-RGB) 
% based on the assessment method developed.

% Consider whether zero-meaning and log-ing are beneficial processes (start
% without and then build in)

clc, clear, close all

%% Load Data
[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_extended',...
    'SSF','SS10',...
    'mel_offset',0,...
    'lum','CIE_10');

%% Compute colorimetry

%First compute colorimetry under a mean illuminant, to work out which
%reflectances to exclude
[~, lsri_avIll] = melcomp_colorimetry(mean(T_SPD,2), T_SRF, T_SSF, T_lum, S_sh);
exclRef = excludeSimilarRefs(lsri_avIll);
T_SRF_reduced = T_SRF(:,~exclRef);

[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF_reduced, T_SSF, T_lum, S_sh);

%% Normalise data for std and set mean to 0

% % Testing standardness of distributions
% hk=10;
% figure,hist(lsri_c(1,:),hk)
% figure,hist(lsri_c(2,:),hk)
% 
% figure,hist(log(lsri_c(1,:)),hk)
% figure,hist(log(lsri_c(2,:)),hk)

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

%% Perform corrections

% Do nothing 'correction' 
% -----------------------
% Baseline measurement
output(:,:,:,1) = lsri(1:2,:,:);


% Grey-world assumption:
% ----------------------
% Here, the average chromaticity of the reflectances under a single 
% illuminant will be taken to be the origin, and everything will be denoted
% in reference to that.

gw = zeros(size(lsri(1:2,:,:)));
for i=1:size(T_SPD,2)
    IllEst = mean(lsri(1:2,:,i),2); %Illuminant estimate
    gw(:,:,i) = lsri(1:2,:,i) - IllEst;
end
output(:,:,:,2) = gw;

% % Figure to check operation
% figure, hold on
% scatter(lsri(1,:,end),lsri(2,:,end),'k')
% scatter(IllEst(1),IllEst(2),'r','filled')
% scatter(gw(1,:,end),gw(2,:,end),'ks')
% plot([0,0],[min(ylim),max(ylim)],'k--')
% plot([min(xlim),max(xlim)],[0,0],'k--')
% scatter(0,0,'rs','filled')


% Bright-is-white assumption:
% ---------------------------
% Here, the surface with the highest luminance value will be taken as being
% a neutral chromaticity, and will be shifted to the origin, with
% everything else denoted in reference to that.

BiW = zeros(size(lsri(1:2,:,:)));
maxLumLoc = zeros(size(T_SPD,2),1);
for i=1:size(T_SPD,2)
    Lum = T_lum'*(T_SRF_reduced.*T_SPD(:,i));
    [~,maxLumLoc(i)] = max(Lum);
    IllEst(:,i) = lsri(1:2,maxLumLoc(i),i); %Illuminant estimate
    BiW(:,:,i) = lsri(1:2,:,i) - IllEst(:,i);
end
output(:,:,:,3) = BiW;

% % Figures to check operation
% figure, hold on
% for i=1:500:size(T_SPD,2)
% scatter(lsri(1,:,i),lsri(2,:,i))
% ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex-1; %resets colour order so that colours correspond
% scatter(BiW(1,:,i),BiW(2,:,i),'s')
% end
% plot([0,0],[min(ylim),max(ylim)],'k--')
% plot([min(xlim),max(xlim)],[0,0],'k--')
% figure, 
% hist(maxLumLoc,100)


% Melanopsin based correction
% ---------------------------

l_cal_range = 0:0.001:1.5;
s_cal_range = -2:0.001:-0.15;

[sf_l,sf_s] = melcomp_6_calcsf(lsri, l_cal_range, s_cal_range); %calculates scaling factors
MC = lsri;
MC(1,:) = MC(1,:)+sf_l*MC(4,:);
MC(2,:) = MC(2,:)+sf_s*MC(4,:);
output(:,:,:,4) = MC(1:2,:,:);

%% Summary figure

pltc_alt = repmat(jet(size(T_SRF_reduced,2))',1,1,size(T_SPD,2));
titles = {'Baseline','Grey World','Bright-is-white','Melanopsin Correction'};

figure,
for i=1:size(output,4)
    subplot(sqrt(size(output,4)),sqrt(size(output,4)),i)
    scatter(reshape(output(1,:,:,i),[],1),reshape(output(2,:,:,i),[],1),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    title(titles{i})
    axis equal
end

%% Judgement time

for i=1:size(output,4)
mark(i) = KMeansMark(squeeze(output(:,:,:,i)));
end

mark




