clear
clc

%%

S = [390 1 441]; 
SToWls(S)

T_quantal_PTB = ComputeCIEConeFundamentals(S,10,32,3);
T_energy_PTB = EnergyToQuanta(S,T_quantal_PTB')';

%%

figure, plot(T_energy_PTB')
%figure, plot(T_quantal_PTB')