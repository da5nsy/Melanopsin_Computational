clear
clc

%% Generate cone fundamentals

S = [390 1 441]; 
%SToWls(S)

T_quantal_PTB = ComputeCIEConeFundamentals(S,10,32,3);
%figure, plot(T_quantal_PTB')

%%

%Scaling factors from Eq. 8.5 from CIE 170-2:2015
sf_q10 = [0.67390486, 0.35656342, 0.06834828]; %quantal 10
sf_e10 = [0.69283932, 0.34967567, 0.05547858]; %energy 10

%% Test by computing EE co-ordinates

EEenergy = ones(length(T_quantal_PTB),1); %equi-energy white
EEquanta = EnergyToQuanta(S,EEenergy);

LMS_EE = sf_q10.*(T_quantal_PTB*EEquanta);

%LMS_EE = sum(T_quantal_PTB').*sf;

l = LMS_EE(1)/(LMS_EE(1)+LMS_EE(2))
m = LMS_EE(2)/(LMS_EE(1)+LMS_EE(2));
s = LMS_EE(3)/(LMS_EE(1)+LMS_EE(2))

% l =
% 
%     0.7091
% 
% 
% s =
% 
%     0.0207


% From CIE 170-2:2015, pg 12:

% "The MacLeod–Boynton chromaticity diagram for 10° field size is shown in Figure 8.3. In
% particular, Illuminant E (the equi-energetic spectrum) is represented at
% (0,699 237; 0,025 841)."

%% Test by computing spectral co-ordinates

l_spectral = sf_q10(1)*T_quantal_PTB(1,:)./(sf_q10(1)*T_quantal_PTB(1,:)+sf_q10(2)*T_quantal_PTB(2,:));
m_spectral = sf_q10(2)*T_quantal_PTB(2,:)./(sf_q10(1)*T_quantal_PTB(1,:)+sf_q10(2)*T_quantal_PTB(2,:));
s_spectral = sf_q10(3)*T_quantal_PTB(3,:)./(sf_q10(1)*T_quantal_PTB(1,:)+sf_q10(2)*T_quantal_PTB(2,:));

figure,
scatter(l_spectral,s_spectral)
ylim([0,1])


