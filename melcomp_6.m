clear
clc

%%

S = [390 1 441]; 
%SToWls(S)

T_quantal_PTB = ComputeCIEConeFundamentals(S,10,32,3);
T_energy_PTB = EnergyToQuanta(S,T_quantal_PTB')';

%figure, plot(T_energy_PTB')

%%

%Scaling factors from Eq. 8.5 from CIE 170-2:2015
sf = [0.69283932, 0.34967567, 0.05547858];

LMS_EE = sum(T_energy_PTB').*sf;

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


