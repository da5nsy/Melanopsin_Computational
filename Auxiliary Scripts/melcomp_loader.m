function [T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(varargin)

%%
 
%clear, clc, close all
%varargin = {'SPD','Granada','SRF','Vrhel_nat_1','SSF','SP','mel_offset',0}; %just for testing

%% Parse inputs
expectedSPD = {'Granada','D-series','Granada_sub'};
expectedSRF = {'Vrhel_nat_1','Vrhel_nat_2','Vrhel_nat_extended','Vrhel_full','Foster'};
expectedSSF = {'SS10','SP'};
expectedlum = {'CIE_10','SP'};

default_SPD = expectedSPD{1};
default_SRF = expectedSRF{1};
default_SSF = expectedSSF{1};
default_mel_offset = 0;
default_lum = expectedlum{1};

p = inputParser;
addParameter(p,'SPD',default_SPD, @(x) any(validatestring(x,expectedSPD)));
addParameter(p,'SRF',default_SRF, @(x) any(validatestring(x,expectedSRF)));
addParameter(p,'SSF',default_SSF, @(x) any(validatestring(x,expectedSSF)));
addParameter(p,'mel_offset',default_mel_offset);
addParameter(p,'lum',default_lum, @(x) any(validatestring(x,expectedlum)));

parse(p,varargin{:});

%% T_SPD

if strcmp(p.Results.SPD,'Granada') || strcmp(p.Results.SPD,'Granada_sub')
    load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
    T_SPD = final; clear final
    S_SPD = [300,5,161];
    
    if strcmp(p.Results.SPD,'Granada_sub')
        T_SPD = T_SPD(:,1:20:end);
    end
    
    % http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
    % From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
    %       spectroradiometric characteristics of narrow-field-of-view
    %       clear skylight in Granada, Spain" (2001)
end

if strcmp(p.Results.SPD,'D-series')
    D_CCT=1./linspace(1/3600,1/25000,20); %non-linear range, aiming to better reproduce observed variation
    load B_cieday.mat B_cieday S_cieday
    T_SPD = GenerateCIEDay(D_CCT,[B_cieday]); %these appear to be linearly upsampled from 10nm intervals (see 'cieday investigation.m' https://github.com/da5nsy/General-Purpose-Functions/blob/3ee587429e9c4f3dd52d64acd95acf82d7e05f47/cieday%20investigation.m)
    T_SPD = T_SPD./max(T_SPD); %normalise
    S_SPD = S_cieday;
end

%% T_SRF

if strcmp(p.Results.SRF,'Vrhel_nat_1')
    load sur_vrhel.mat sur_vrhel S_vrhel   
    refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
    T_SRF=sur_vrhel(:,refs);
    S_SRF = S_vrhel;
end

if strcmp(p.Results.SRF,'Vrhel_nat_2')
    load sur_vrhel.mat sur_vrhel S_vrhel   
    refs=[38, 15, 134, 137, 138, 65, 19, 24, 140, 26];
    T_SRF=sur_vrhel(:,refs);
    S_SRF = S_vrhel;
end

if strcmp(p.Results.SRF,'Vrhel_nat_extended')
    load sur_vrhel.mat sur_vrhel S_vrhel   
    refs=[1:44,65,69,81:154];
    T_SRF=sur_vrhel(:,refs);
    S_SRF = S_vrhel;
end

if strcmp(p.Results.SRF,'Vrhel_full')
    load sur_vrhel.mat sur_vrhel S_vrhel   
    T_SRF=sur_vrhel;
    S_SRF = S_vrhel;
end

if strcmp(p.Results.SRF,'Foster') %this is a little messy. Need to check whether this refers to uncropped/cropped images !!!!!!!!!!
    base = 'C:\Users\cege-user\Documents\Large data\Foster Images\';
    %     for i=1:4 %2002 images
    %         ims(i)=load([base, '2002\scene',num2str(i),'.mat']); %imageS
    %     end
    %     %2004 images
    ims(5)=load([base,'2004\scene1\ref_crown3bb_reg1.mat']);
    ims(6)=load([base,'2004\scene2\ref_ruivaes1bb_reg1.mat']);
    %     ims(7)=load([base,'2004\scene3\ref_mosteiro4bb_reg1.mat']);
    %     ims(8)=load([base,'2004\scene4\ref_cyflower1bb_reg1.mat']);
    %     ims(9)=load([base,'2004\scene5\ref_cbrufefields1bb_reg1.mat']);
    
    ims(6).reflectances = ims(6).reflectances(:,1:end-100,:);
    [r, c, w] = size(ims(6).reflectances);
    T_SRF = reshape(ims(6).reflectances, r*c, w);
    %     for i=6%length(ims) %add this loop to append multiple images
    %         [r, c, w] = size(ims(i).reflectances);
    %         T_refs = [T_refs; reshape(ims(i).reflectances, r*c, w)];
    %     end
    %S_refs=[410,10,31]; %for 2002 data
    S_SRF=[400,10,33];%for 2004 data (except #9)
end

%% T_SSF

if or(...
        and(strcmp(p.Results.SSF,'SS10'),strcmp(p.Results.lum,'SP')),...
        and(strcmp(p.Results.SSF,'SP'),  strcmp(p.Results.lum,'CIE_10')))
    error('Mismatch between SSF and lum')
end

if strcmp(p.Results.SSF,'SS10')
    % Stockman & Sharpe 10deg, for use with MacLeod Boynton diagram, 
    % as per CIE 170-2:2015
    load T_cones_ss10.mat T_cones_ss10 S_cones_ss10
    T_SSF = T_cones_ss10;
    S_SSF = S_cones_ss10;
end

if strcmp(p.Results.lum,'CIE_10')
    load T_CIE_Y10.mat T_CIE_Y10 S_CIE_Y10
    T_lum = T_CIE_Y10;
    S_lum = S_CIE_Y10;
end

if strcmp(p.Results.SSF,'SP')
    % Smith-Pokorny, for use with original type MacLeod Boynton diagram
    load T_cones_sp.mat T_cones_sp S_cones_sp
    T_SSF = T_cones_sp;
    S_SSF = S_cones_sp;
end

if strcmp(p.Results.lum,'SP')
    T_lum = [0.6373 0.3924 0]*T_cones_sp; %figures from PTB LMSToMacBoyn.m
    S_lum = S_cones_sp;
end

load T_rods T_rods S_rods
load T_melanopsin T_melanopsin S_melanopsin
S_melanopsin_big = [200,1,800]; %Increasing the range of S_melanopsin to allow for shifting without affecting S_sh
T_melanopsin_big = SplineCmf(S_melanopsin,T_melanopsin,S_melanopsin_big);

%% Mel offset
if p.Results.mel_offset
    if p.Results.mel_offset < -S_melanopsin_big(1)-1
        error('You have gone too low on the melanopsin')
    else
        S_melanopsin_big(1) = S_melanopsin_big(1)+ p.Results.mel_offset;
    end
end

%% Reduce all data down to the common range/interval

S_sh = [max([S_SPD(1),S_SRF(1),S_SSF(1),S_rods(1),S_melanopsin_big(1)]),...
    max([S_SPD(2),S_SRF(2),S_SSF(2),S_rods(2),S_melanopsin_big(2)]),...
    min([S_SPD(3),S_SRF(3),S_SSF(3),S_rods(3),S_melanopsin_big(3)])];
%S_shared: work out what the lowest common denominator for the range/interval of the data is

T_SPD  = SplineSpd(S_SPD,T_SPD,S_sh);
T_SRF  = SplineSrf(S_SRF,T_SRF,S_sh,1); %ended with same value
T_SSF  = SplineCmf(S_SSF,T_SSF,S_sh,1)';
T_rods = SplineCmf(S_rods,T_rods,S_sh,0)';
T_mel  = SplineCmf(S_melanopsin_big,T_melanopsin_big,S_sh,0)';
T_lum  = SplineCmf(S_lum,T_lum,S_sh,0)';

% combine sensitivity vectors
T_SSF=[T_SSF,T_rods,T_mel];

end