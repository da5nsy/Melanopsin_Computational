% Practical version of the abstract formulae in 

% Maloney, L. T. Computational Approaches to Color Constancy. (Stanford, 1984).
% and 
% SURFACE COLOUR PERCEPTION AND ENVIRONMENTAL CONSTRAINTS, laurence t. maloney
% in Mausfeld, R. & Heyer, D. Colour Perception: Mind and the Physical World. 
% (Oxford University Press, USA, 2004).

clear, clc, close all

%% Load basis functions and define
% 1. Basis functions
% 2. Basis weights
% 3. Full spectra

load B_cieday
load B_vrhel

[E,~,~,e(2),e(3)] = GenerateCIEDay(6500,B_cieday); 
e(1) = 1;
e=e';
%Calculates an SPD, and weights of the basis functions also the weights of
%the basis functions needed to create that (E(spd) = B_cieday*e)
% check = isequal(E, B_cieday*e);

B_vrhel = B_vrhel(:,1:3); %Let's pretend it only has 3 dimensions, for simplicity
s = [-1;0.2;0.6]; %random picks
S = B_vrhel*s;

figure, hold on
plot(SToWls(S_cieday),E/max(E),'b','LineWidth',4)
for i =1:3
plot(SToWls(S_cieday),B_cieday(:,i)*e(i)/max(E),'b--','LineWidth',4-i)
end

plot(SToWls(S_vrhel),S/max(S),'r','LineWidth',5)
for i =1:3
plot(SToWls(S_vrhel),B_vrhel(:,i)*s(i)/max(S),'r--','LineWidth',4-i)
end

plot([min(xlim),max(xlim)],[0,0],'k:')
legend({'E','B_cieday(1) * e(1)','B_cieday(2) * e(2)','B_cieday(3) * e(3)'...
    'S','B_vrhel(1) * s(1)','B_vrhel(2) * s(2)','B_vrhel(3) * s(3)'},'Location','best','Interpreter', 'none')


%% 


load T_cones_sp

%%

% (eq numbers from second ref)

% Eq 9.3

L = E*S;
%p = L*R;

