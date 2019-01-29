function  pc = melcomp_6_looper(offset,norm,plt) 

if ~exist('offset','var') %do this properly with nargin !!!!!!!!!!!!!
    disp('No offset passed. Setting offset to 0.')
    offset = 0;
    disp('No norm command passed. Setting norm to 1 (positive)')
    offset = 1;
end

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

S_mel(1)=S_mel(1)+offset; % __________THIS IS THE IMPORTANT BIT___________

% Load reflectance data

load sur_vrhel
refs=[38, 15, 134, 137, 138, 65, 19, 24, 140, 26];
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
T_SPD = T_SPD(:,1:20:end);
S_SPD=[300,5,161];

% reduce all data down to the common range/interval
S_sh = [max([S_SPD(1),S_SRF(1),S_SSF(1),S_rods(1),S_mel(1)]),...
    max([S_SPD(2),S_SRF(2),S_SSF(2),S_rods(2),S_mel(2)]),...
    min([S_SPD(3),S_SRF(3),S_SSF(3),S_rods(3),S_mel(3)])]; 
%S_shared: work out what the lowest common denominator for the range/interval of the data is

T_SPD = SplineSpd(S_SPD,T_SPD,S_sh);
T_SRF = SplineSrf(S_SRF,T_SRF,S_sh,1); %ended with same value
T_SSF  = SplineCmf(S_SSF,T_SSF,S_sh)';
T_rods = SplineCmf(S_rods,T_rods,S_sh)';
T_mel  = SplineCmf(S_mel,T_mel,S_sh)';
[S_SPD, S_SRF, S_SSF, S_rods, S_mel] = deal(S_sh);

% combine sensitivity vectors
T_LMSRI=[T_SSF,T_rods,T_mel];
S_LMSRI=S_sh;

%

sf_10 = [0.69283932, 0.34967567, 0.05547858]; %energy 10deg from CIE 170-2:2015
T_lum = sf_10(1)*T_SSF(:,1)+sf_10(2)*T_SSF(:,2);

%

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

%% correct axes for variance

% figure,
% scatter3(lsri(1,:),lsri(2,:),lsri(4,:));

lsri = lsri(:,:);

if norm
    for i=1:4
        lsri(i,:) = (lsri(i,:) - mean(lsri(i,:)))./std(lsri(i,:));
    end
end


%% pca

lsi = lsri([1,2,4],:);

[pc.coeff, pc.score, pc.latent, pc.tsquared, pc.explained, pc.mu] = pca(lsi');

%%

if plt
    figure,
    scatter3(lsri(1,:),lsri(2,:),lsri(4,:));
end
    

end