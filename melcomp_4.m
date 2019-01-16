clc, clear, close all

%% SPDs (Spectral Power Distributions)

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
T_SPD=final; clear final
S_SPD=[300,5,161];

T_SPD = T_SPD(:,1:100:end);

%% SRFs (Spectral Reflectance Functions)

load sur_vrhel
T_SRF = sur_vrhel(:,[1:44,65,69,81:154]); %natural only
S_SRF = S_vrhel;
clear sur_vrhel S_vrhel

%% SSF (Spectral Sensitivity Functions)

load T_cones_sp T_cones_sp S_cones_sp
T_SSF = T_cones_sp;
S_SSF = S_cones_sp;
clear T_cones_sp S_cones_sp

load T_rods T_rods S_rods
load T_melanopsin T_melanopsin S_melanopsin
T_mel = SplineCmf(S_melanopsin, T_melanopsin, S_melanopsin - [10, 0, 0],1); %Increasing the range of this function in case it ends up limiting the range of S_sh, and shorten variable names
S_mel = S_melanopsin - [10, 0, 0];
clear S_melanopsin T_melanopsin

%% Interpolate functions to match range and interval

%reduce all data down to the common range/interval
S_sh = [max([S_SPD(1),S_SRF(1),S_SSF(1)]),max([S_SPD(2),S_SRF(2),S_SSF(2)]),min([S_SPD(3),S_SRF(3),S_SSF(3)])]; %S_shared: work out what the lowest common denominator for the range/interval of the data is

T_SPD = SplineSpd(S_SPD,T_SPD,S_sh,1); % extend == 1: Cubic spline, extends with last value in that direction
T_SRF = SplineSrf(S_SRF,T_SRF,S_sh,1);
T_SSF  = SplineCmf(S_SSF,T_SSF,S_sh,1)';
T_rods = SplineCmf(S_rods,T_rods,S_sh)'; %extend? !!!!!!!!!!
T_mel  = SplineCmf(S_mel,T_mel,S_sh)'; %extend? !!!!!!!!!!
[S_SPD, S_SRF, S_SSF, S_rods, S_mel] = deal(S_sh);

% combine sensitivity vectors
T_LMSRI=[T_SSF,T_rods,T_mel];
S_LMSRI=S_sh;

%% Calculate signals

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_SRF.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_LMSRI'*T_rad(:,:,i);
    
    LMSRI_ill(:,:,i) = T_LMSRI'*T_SPD(:,i); %purposefully left the second dim in for now, as a reminder that it essentially represents a single reflectance
end

LMSRI_ill = repmat(LMSRI_ill,1,120);

slm_obj = LMSRI(3,:,:)./(LMSRI(1,:,:));
si_obj = LMSRI(3,:,:)./LMSRI(5,:,:);

slm_ill = LMSRI_ill(3,:,:)./(LMSRI_ill(1,:,:)+LMSRI_ill(2,:,:));

%%

figure,
subplot(1,2,1)
scatter(slm_obj(:),slm_ill(:),'k.')
xlim([0 2])
xlabel('slm_obj','Interpreter','None')
ylabel('slm_ill','Interpreter','None')

subplot(1,2,2)
scatter(si_obj(:),slm_ill(:),'k.')
xlim([0 2])
xlabel('si_obj','Interpreter','None')

%%

figure,
subplot(1,2,1)
scatter(mean(slm_obj,2),mean(slm_ill,2),'k.')
xlim([0 2])
xlabel('slm_obj','Interpreter','None')
ylabel('slm_ill','Interpreter','None')

subplot(1,2,2)
scatter(mean(si_obj,2),mean(slm_ill,2),'k.')
xlim([0 2])
xlabel('si_obj','Interpreter','None')