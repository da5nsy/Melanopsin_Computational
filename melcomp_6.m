clear
clc

%% Generate cone fundamentals

load T_cones_ss10.mat
figure, plot(SToWls(S_cones_ss10),T_cones_ss10')

%%

%Scaling factors from Eq. 8.5 from CIE 170-2:2015

sf_e10 = [0.69283932, 0.34967567, 0.05547858]; %energy 10

%% Test by computing EE co-ordinates

EEenergy = ones(length(T_cones_ss10),1); %equi-energy white

LMS_EE_e = sf_e10.*(T_cones_ss10*EEenergy);

l_e = LMS_EE_e(1)/(LMS_EE_e(1)+LMS_EE_e(2));
m_e = LMS_EE_e(2)/(LMS_EE_e(1)+LMS_EE_e(2));
s_e = LMS_EE_e(3)/(LMS_EE_e(1)+LMS_EE_e(2));

% From CIE 170-2:2015, pg 12:

% "The MacLeod–Boynton chromaticity diagram for 10° field size is shown in Figure 8.3. In
% particular, Illuminant E (the equi-energetic spectrum) is represented at
% (0,699 237; 0,025 841)."

%% Test by computing spectral co-ordinates

l_spectral_e = sf_e10(1)*T_cones_ss10(1,:)./(sf_e10(1)*T_cones_ss10(1,:)+sf_e10(2)*T_cones_ss10(2,:));
m_spectral_e = sf_e10(2)*T_cones_ss10(2,:)./(sf_e10(1)*T_cones_ss10(1,:)+sf_e10(2)*T_cones_ss10(2,:));
s_spectral_e = sf_e10(3)*T_cones_ss10(3,:)./(sf_e10(1)*T_cones_ss10(1,:)+sf_e10(2)*T_cones_ss10(2,:));

figure,
scatter(l_spectral_e,s_spectral_e)
ylim([0,1])


