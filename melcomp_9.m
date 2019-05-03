% The goal here is to compare a melanopsin based correction with other
% basic corrections (do nothing, grey world, bright-is-white/max-RGB)
% based on the assessment method developed.

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

% Calculate luminance for BiW
Lum = zeros(size(T_SRF_reduced,2),size(T_SPD,2));
for i=1:size(lsri,3)
    Lum(:,i) = T_lum'*(T_SRF_reduced.*T_SPD(:,i));
end

rng(1)
pcSurfRange = 5:5:100; %percent surfaces range
nMethods = 4;
mark = zeros(nMethods,length(pcSurfRange));
for pcSurf = 1:length(pcSurfRange)
    nSurf = round((pcSurfRange(pcSurf)/100) * size(T_SRF_reduced,2));
    
    lsri2 = zeros(size(lsri,1),nSurf,size(T_SPD,2));
    Lum2 = zeros(nSurf,size(T_SPD,2));
    sel = zeros(nSurf,size(T_SPD,2));
    for i = 1:size(T_SPD,2)
        sel(:,i) = randperm(size(T_SRF_reduced,2),nSurf); %selection
        lsri2(:,:,i) = lsri(:,sel(:,i),i);
        Lum2(:,i) = Lum(sel(:,i),i); %!!!!!!!!!!!!!!!!
    end
    
    % Perform corrections
    output = performCC(lsri2,Lum2);
    
    % Score corrections
    for i=1:nMethods
        mark(i,pcSurf) = KMeansMark(squeeze(output(:,:,:,i)),size(T_SRF_reduced,2),sel);
    end
    
    disp(pcSurfRange(pcSurf))
    
%     plot(pcSurfRange,mark','.')
%     xlim([0 100])
%     ylim([0 1])
%     drawnow
end


%% Summary figure

pltc_alt = repmat(jet(size(T_SRF_reduced,2))',1,1,size(T_SPD,2));
titles = {'Baseline','Grey World','Bright-is-white','Melanopsin Correction'};


figure, hold on
plot(pcSurfRange,mark')
legend(titles,'Location','best')
plot(pcSurfRange,1./round((pcSurfRange/100) * size(T_SRF_reduced,2)),'k','DisplayName','Chance')
xlim([0 100])
ylim([0 1])

figure,
for i=1:size(output,4)
    subplot(sqrt(size(output,4)),sqrt(size(output,4)),i)
    scatter(reshape(output(1,:,:,i),[],1),reshape(output(2,:,:,i),[],1),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    title(titles{i})
    axis equal
end



