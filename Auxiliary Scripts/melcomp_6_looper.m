function  pc = melcomp_6_looper(varargin) 

% clear, clc, close all
% varargin = {'mel_offset',0,'norm',1,'plt',0}; %just for testing

default_mel_offset = 0;
default_norm       = 1;
default_plt        = 0;

expectedSPD = {'Granada_sub','Granada','D-series'};
expectedSRF = {'Vrhel_nat_1','Vrhel_nat_2','Vrhel_full','Foster'};
expectedSSF = {'SS10','SP'};
expectedlum = {'CIE_10','SP'};
default_SPD = expectedSPD{1};
default_SRF = expectedSRF{2}; %Note - this is different behaviour to melcomp_loader since melcomp_6 uses Vrhel_nat_2 by default
default_SSF = expectedSSF{1};
default_lum = expectedlum{1};

p = inputParser;
addParameter(p,'mel_offset',default_mel_offset);
addParameter(p,'norm',default_norm);
addParameter(p,'plt',default_plt);
addParameter(p,'LMSToMacBoyn_old',0);
addParameter(p,'SPD',default_SPD, @(x) any(validatestring(x,expectedSPD)));
addParameter(p,'SRF',default_SRF, @(x) any(validatestring(x,expectedSRF)));
addParameter(p,'SSF',default_SSF, @(x) any(validatestring(x,expectedSSF)));
addParameter(p,'lum',default_lum, @(x) any(validatestring(x,expectedlum)));

parse(p,varargin{:});

%% Data

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD',p.Results.SPD,...
    'SRF',p.Results.SRF,...
    'SSF',p.Results.SSF,...
    'mel_offset',p.Results.mel_offset,...
    'lum',p.Results.lum);

%

T_rad = zeros([S_sh(3),size(T_SRF,2),size(T_SPD,2)]);
LMSRI = zeros([size(T_SSF,2),size(T_SRF,2),size(T_SPD,2)]);
lsri  = zeros([4,size(T_SRF,2),size(T_SPD,2)]);
t_r   = zeros([2,size(T_SRF,2),size(T_SPD,2)]); %t for temp
t_i   = zeros([2,size(T_SRF,2),size(T_SPD,2)]); %t for temp

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_SRF.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_SSF'*T_rad(:,:,i);
    if p.Results.LMSToMacBoyn_old
        lsri(1:2,:,i) = LMSToMacBoyn_old(LMSRI(1:3,:,i),T_SSF(:,1:3)',T_lum');
        t_r(:,:,i)    = LMSToMacBoyn_old(LMSRI([1,2,4],:,i),[T_SSF(:,1:2)';T_SSF(:,4)'],T_lum');
        t_i(:,:,i)    = LMSToMacBoyn_old(LMSRI([1,2,5],:,i),[T_SSF(:,1:2)';T_SSF(:,5)'],T_lum');
    else
        lsri(1:2,:,i) = LMSToMacBoyn(LMSRI(1:3,:,i),T_SSF(:,1:3)',T_lum');
        t_r(:,:,i)    = LMSToMacBoyn(LMSRI([1,2,4],:,i),[T_SSF(:,1:2)';T_SSF(:,4)'],T_lum');
        t_i(:,:,i)    = LMSToMacBoyn(LMSRI([1,2,5],:,i),[T_SSF(:,1:2)';T_SSF(:,5)'],T_lum');
    end
end
lsri(3,:,:) = t_r(2,:,:); clear t_r
lsri(4,:,:) = t_i(2,:,:); clear t_i

%% correct axes for variance

% figure,
% scatter3(lsri(1,:),lsri(2,:),lsri(4,:));

lsri = lsri(:,:);

if p.Results.norm
    for i=1:4
        lsri(i,:) = (lsri(i,:) - mean(lsri(i,:)))./std(lsri(i,:));
    end
end


%% pca

lsi = lsri([1,2,4],:);

[pc.coeff, pc.score, pc.latent, pc.tsquared, pc.explained, pc.mu] = pca(lsi');

%%

if p.Results.plt
    %figure,
    scatter3(lsri(1,:),lsri(2,:),lsri(4,:),'k','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6);
end
    

end