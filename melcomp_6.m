% Testing the calculation of MacLeod-Boynton chromaticities based on CIE
% 170-1:2006 and 170-2:2015.

%% Generate cone fundamentals

load T_cones_ss10.mat % Stockman & Sharpe cone fundamentals (10deg)
% % From PsychToolbox. Also available from CIE 170-1:2006, table 6.2, or cvrl.org

figure, plot(SToWls(S_cones_ss10),T_cones_ss10')

%%

% % Scaling factors from Eq. 8.5 from CIE 170-2:2015, pg 8

sf_10 = [0.69283932, 0.34967567, 0.05547858]; %energy 10deg

%% Test by computing EE co-ordinates

EE = ones(length(T_cones_ss10),1); % define equi-energy (EE) white

LMS_EE = sf_10.*(T_cones_ss10*EE)'; %calcaulte EE tristimulus values

l_EE = LMS_EE(1)/(LMS_EE(1)+LMS_EE(2)); %calculate EE chromaticity co-ordinates
m_EE = LMS_EE(2)/(LMS_EE(1)+LMS_EE(2)); %unused
s_EE = LMS_EE(3)/(LMS_EE(1)+LMS_EE(2));

% % From CIE 170-2:2015, pg 12:
% % 
% % "The MacLeod–Boynton chromaticity diagram for 10° field size is shown in Figure 8.3. In
% % particular, Illuminant E (the equi-energetic spectrum) is represented at
% % (0,699 237; 0,025 841)."

% disp(l_EE)
% disp(s_EE)

%% Test by computing spectral co-ordinates

% calculate chromaticity co-ordinates of spectral locus
l_spectral = sf_10(1)*T_cones_ss10(1,:)./(sf_10(1)*T_cones_ss10(1,:)+sf_10(2)*T_cones_ss10(2,:));
m_spectral = sf_10(2)*T_cones_ss10(2,:)./(sf_10(1)*T_cones_ss10(1,:)+sf_10(2)*T_cones_ss10(2,:)); %unused
s_spectral = sf_10(3)*T_cones_ss10(3,:)./(sf_10(1)*T_cones_ss10(1,:)+sf_10(2)*T_cones_ss10(2,:));

figure, hold on
scatter(l_spectral,s_spectral,'.') %spectral locus
xlim([0.4 1])
ylim([0,1])

scatter(l_EE,s_EE,'*') %EE white chromaticity

