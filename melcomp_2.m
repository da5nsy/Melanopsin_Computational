function [pc, LMSRI] = melcomp_2(varargin)


%% Pre-flight checks

try
    nargin;
catch
    clear, clc, close all
    varargin = {'SPD','Granada_sub','SRF','Vrhel_nat_2'};
end

default_mel_offset = 0;
default_Z_ax       = 9;
default_plt        = '';

expectedSPD = {'Granada_sub','Granada','D-series'};
expectedSRF = {'Vrhel_nat_1','Vrhel_nat_2','Vrhel_full','Foster'};
expectedSSF = {'SS10','SP'};
expectedlum = {'CIE_10','SP'};
default_SPD = expectedSPD{3};
default_SRF = expectedSRF{1};
default_SSF = expectedSSF{2};
default_lum = expectedlum{2};

p = inputParser;
addParameter(p,'mel_offset',default_mel_offset);
addParameter(p,'SPD',default_SPD, @(x) any(validatestring(x,expectedSPD)));
addParameter(p,'SRF',default_SRF, @(x) any(validatestring(x,expectedSRF)));
addParameter(p,'SSF',default_SSF, @(x) any(validatestring(x,expectedSSF)));
addParameter(p,'lum',default_lum, @(x) any(validatestring(x,expectedlum)));

addParameter(p,'Z_ax',default_Z_ax);
addParameter(p,'plt',default_plt);

parse(p,varargin{:});

%% Load Data

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD',p.Results.SPD,...
    'SRF',p.Results.SRF,...
    'SSF',p.Results.SSF,...
    'mel_offset',p.Results.mel_offset,...
    'lum',p.Results.lum);

%% Compute colorimetry

[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

%compute colours for display
pltc_alt = repmat(jet(size(T_SRF,2))',1,1,size(T_SPD,2));
rng(7); pltc_alt=pltc_alt(:,randperm(size(T_SRF,2)),:); %best combo for differentiating close chromaticities

%% Plot third dimension

plt_3D = 0; % hard-code, 0 = off, 1 = on

plt_locus = 0; % plot spectral locus in the MB diagram, 0 = off, 1 = on

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

if plt_3D || strcmp(p.Results.plt,'3D')
    
    %figure,
   
    if ismember(p.Results.Z_ax,1:5)
        t_Z = LMSRI(p.Results.Z_ax,:,:); %temp Z
    elseif ismember(p.Results.Z_ax,6:9)
        t_Z = lsri(p.Results.Z_ax-5,:,:);
    elseif p.Results.Z_ax == 10
        t_Z = LMSRI(1,:,:)+LMSRI(2,:,:);
    elseif p.Results.Z_ax == 11
        t_Z = 0.6373*LMSRI(1,:,:)+0.3924*LMSRI(2,:,:);
    elseif p.Results.Z_ax == 12
        t_Z = lsri(3,:,:)+lsri(4,:,:);
    end    
    
    scatter3(lsri(1,:),lsri(2,:),t_Z(1,:),[],pltc_alt(:,:)','v','filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)
    hold on
    zlabel(plt_lbls{p.Results.Z_ax})

    if plt_locus
        MB_locus = LMSToMacBoyn(T_SSF(:,1:3)',T_SSF(:,1:3)',T_lum');
        %plot(MB_locus(1,:),MB_locus(2,:))
        fill([MB_locus(1,5:65),MB_locus(1,5)],[MB_locus(2,5:65),MB_locus(2,5)],'k','LineStyle','none','FaceAlpha','0.1')
    end
    
    grid on
    %axis equal
    %xlim([0 1]), ylim([0 1])
    xlabel('l'),ylabel('s'), 
    %title(plt_lbls{p.Results.Z_ax})
    %view(3) 
    view(188,46)
end

%% Correction through rotation

plt_CTR = 0;

%rotation matrix
ang1  = 0.27; %angle in radians, just eyeballed (for Granada data)
ang2  = -0.90;

rm=...
    [cos(ang2),0,0,sin(ang2);...
    0,cos(ang1),0,sin(ang1);...
    0,0,1,0;...
    -sin(ang2),-sin(ang1),0,cos(ang1)+cos(ang2)]; 

% %Single dimension, only calibrate s
% rm=...
%     [1,0,0,0;...
%     0,cos(ang1),0,sin(ang1);...
%     0,0,1,0;...
%     0,-sin(ang1),0,cos(ang1)]; 

%apply rotation
lsri_rs=lsri(:,:)'*rm;

if plt_CTR || strcmp(p.Results.plt,'CTR')
    figure,
    scatter3(lsri(1,:),lsri(2,:),lsri(4,:),[],pltc_alt(:,:)','v','filled')
    hold on
    grid on
    % %xlim([0 1]), ylim([-1 1]), zlim([0,2])
    scatter3(lsri_rs(:,1),lsri_rs(:,2),lsri_rs(:,4),[],pltc_alt(:,:)','^','filled')
    
    legend({'Original','Rotated'},'Location','best')
    xlabel('l'),ylabel('s'),zlabel('i');
    
    view(0,90)
    %view(90,0)
end

%% Correction through subtractive shift

plt_CTSS= 0;

lsri_ss = lsri; %shifted

lsri_ss(1,:) = lsri(1,:)-(lsri(4,:)*-1.2);
lsri_ss(2,:) = lsri(2,:)-(lsri(4,:)*0.27);

if plt_CTSS
    figure, hold on, 
    %axis equal, 
    grid on
    scatter3(lsri(1,:),lsri(2,:),lsri(4,:),[],pltc_alt(:,:)','v','filled')
    scatter3(lsri_ss(1,:),lsri_ss(2,:),lsri_ss(4,:),[],pltc_alt(:,:)','^','filled')
    
    legend({'Original','Shifted'},'Location','best')
    xlabel('l'),ylabel('s2'),zlabel('i2');
    
    %view(90,0)
    ylim([-0.02 0.04])
end

%% Correction through multiplicative shift

% WORK IN PROGRESS %
% This runs the risk of just making the values really tiny, aka Ford Model
% Colour Constancy.
% I suspect that this is not a feasible method.

plt_CTMS = 0;

lsri_ms = lsri; %shifted

%lsri_ms(2,:) = lsri(2,:).*(1./(lsri(4,:)/max(lsri(4,:))));
%lsri_ms(2,:) = lsri(2,:).*((max(lsri(4,:))-lsri(4,:)));

lsri_ms(2,:) = lsri(2,:).*lsri(4,:)*10;

if plt_CTMS
    figure, hold on, 
    %axis equal, 
    grid on
    scatter3(lsri(1,:),lsri(2,:),lsri(4,:),[],pltc_alt(:,:)','v','filled')
    scatter3(lsri_ms(1,:),lsri_ms(2,:),lsri_ms(4,:),[],pltc_alt(:,:)','^','filled')
    
    legend({'Original','Shifted'},'Location','best')
    xlabel('l'),ylabel('s2'),zlabel('i');
    
    %view(90,0)
end

%% PCA of signals

lsi = lsri([1,2,4],:);

[pc.coeff, pc.score, pc.latent, pc.tsquared, pc.explained, pc.mu] = pca(lsi');

end