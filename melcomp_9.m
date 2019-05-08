function [mark] = melcomp_9(varargin)

% The goal here is to compare a melanopsin based correction with other
% basic corrections (do nothing, grey world, bright-is-white/max-RGB)
% based on the assessment method developed.

if ~exist('varargin','var'), clear, clc, close all, varargin = {}; end
%If we're running this as a script rather than a function, this line sets
% a fresh slate and makes sure we use the defaults below (it acts as if
% we've called it as a function and not passed anything to it)

default_plt = 'all'; %plotting requests
default_pcSurfRange = 12:7:100; %percent surfaces range
default_mel_offset = 0;

p = inputParser;
addParameter(p,'plt',default_plt);
addParameter(p,'pcSurfRange',default_pcSurfRange);
addParameter(p,'mel_offset',default_mel_offset);

parse(p,varargin{:});

%% Load Data
[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_extended',...
    'SSF','SS10',...
    'mel_offset',p.Results.mel_offset,...
    'lum','CIE_10');

%% Compute colorimetry

%First compute colorimetry under a mean illuminant, to work out which
%reflectances to exclude
[~, lsri_avIll] = melcomp_colorimetry(mean(T_SPD,2), T_SRF, T_SSF, T_lum, S_sh);
% Then excludes anything closer than X to another point
exclRef = excludeSimilarRefs(lsri_avIll,0.003);
T_SRF_reduced = T_SRF(:,~exclRef);

[~, lsri] = melcomp_colorimetry(T_SPD, T_SRF_reduced, T_SSF, T_lum, S_sh);

%% Normalise data for std and set mean to 0

% Testing standardness of distributions
% hk=50;
% figure,
% for i=[1,2]
%     subplot(3,2,i*2-1)
%     hist(lsri(i,:),hk), cleanTicks
%     subplot(3,2,i*2)
%     hist(log(lsri(i,:)),hk), cleanTicks
% end
% subplot(3,2,3*2-1)
% hist(lsri(4,:),hk), cleanTicks
% subplot(3,2,3*2)
% hist(log(lsri(4,:)),hk), cleanTicks

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
nMethods = 4;
mark = zeros(nMethods,length(p.Results.pcSurfRange));
output = zeros(2,size(T_SRF_reduced,2),size(T_SPD,2),nMethods,length(p.Results.pcSurfRange));
sel_store = zeros(length(p.Results.pcSurfRange),round((p.Results.pcSurfRange(end)/100) * size(T_SRF_reduced,2)),size(T_SPD,2));
for pcSurf = 1:length(p.Results.pcSurfRange)
    nSurf(pcSurf) = round((p.Results.pcSurfRange(pcSurf)/100) * size(T_SRF_reduced,2));
    
    lsri2 = zeros(size(lsri,1),nSurf(pcSurf),size(T_SPD,2));
    Lum2 = zeros(nSurf(pcSurf),size(T_SPD,2));
    sel = zeros(nSurf(pcSurf),size(T_SPD,2));
    for i = 1:size(T_SPD,2)
        sel(:,i) = randperm(size(T_SRF_reduced,2),nSurf(pcSurf)); %selection
        lsri2(:,:,i) = lsri(:,sel(:,i),i);
        Lum2(:,i) = Lum(sel(:,i),i); %!!!!!!!!!!!!!!!!
    end
    sel_store(pcSurf,1:nSurf(pcSurf),:) = sel;
    
    % Perform corrections
    output(:,1:nSurf(pcSurf),:,:,pcSurf) = performCC(lsri2,Lum2);
    
    % Score corrections
    for i=1:nMethods
        rng(1)
        mark(i,pcSurf) = KMeansMark(squeeze(output(:,1:nSurf(pcSurf),:,i,pcSurf)),size(T_SRF_reduced,2),sel);
    end
    
    disp(p.Results.pcSurfRange(pcSurf))
    
    %     plot(p.Results.pcSurfRange,mark','.')
    %     xlim([0 100])
    %     ylim([0 1])
    %     drawnow
end


%% Summary figure

if or(strcmp(p.Results.plt,'range'),strcmp(p.Results.plt,'all'))
    titles = {'Do nothing','Grey World','Bright-is-White','Melanopsin'};
    markers = {'s:','d:','^:','o:'};
    mfc = hsv(nMethods);
    lw = 3;
    
    figure, hold on
    for i=1:4
        plot(p.Results.pcSurfRange,mark(i,:),markers{i},'Color',mfc(i,:),'MarkerFaceColor',mfc(i,:),'linewidth',lw)
    end
    legend(titles,'Location','Northwest')
    plot(p.Results.pcSurfRange,1./round((p.Results.pcSurfRange/100) * size(T_SRF_reduced,2)),'k','DisplayName','Chance','linewidth',lw-1)
    xlim([0 100])
    ylim([0 1])
    
    xlabel(sprintf('Percentage of surfaces used in each run (/%d)',size(T_SRF_reduced,2)))
    ylabel('K-means-mark (/1)')
    grid on
end


%% Figures for a range of percentages

titles_short = {'DN','GW','BW','MC'};

if or(strcmp(p.Results.plt,'MB'),strcmp(p.Results.plt,'all'))
    jlist = [1,5,length(p.Results.pcSurfRange)];
    figure,
    pltn = 1;
    for i=1:nMethods
        for j = 1:length(jlist)
            subplot(nMethods,length(jlist),pltn)
            gscatter(reshape(output(1,1:nSurf(jlist(j)),:,i,jlist(j)),[],1),...
                reshape(output(2,1:nSurf(jlist(j)),:,i,jlist(j)),[],1),...
                reshape(sel_store(jlist(j),1:nSurf(jlist(j)),:),[],1),...
                [],[],5)
            axis([-5 5 -3 3])
            if jlist(j)~=1
                yticks('')
            else
                ylabel(titles_short{i})
            end
            if i == 1
                title(strcat(string(p.Results.pcSurfRange(jlist(j))),'%(',string(nSurf(end)),')=',string(nSurf(jlist(j)))))
            end
            if i ~= nMethods
                xticks('')
            end
            pltn = pltn+1;
            legend off
        end
    end
end


%% Single figure

if strcmp(p.Results.plt,'single')
    gscatter(reshape(output(1,1:nSurf(length(nSurf)),:,4,length(p.Results.pcSurfRange)),[],1),...
             reshape(output(2,1:nSurf(length(nSurf)),:,4,length(p.Results.pcSurfRange)),[],1),...
             reshape(sel_store(length(nSurf),:,:),[],1),[],[],5)
         legend off
end



end
