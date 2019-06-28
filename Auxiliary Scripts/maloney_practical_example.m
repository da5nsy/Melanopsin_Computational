% Practical version of the abstract formulae in 

% Maloney, L. T., Computational Approaches to Color Constancy. (Stanford, 1984).
%
% and 
%
% Maloney, L. T., SURFACE COLOUR PERCEPTION AND ENVIRONMENTAL CONSTRAINTS
% in Mausfeld, R. & Heyer, D. Colour Perception: Mind and the Physical World. 
% (Oxford University Press, USA, 2004).

% This essentially shows that from knowing the matrix g (which bundles 
% together sensor spectral sensitivy, and bases for illuminant and 
% surfaces) and the weights e (illuminant basis weights), it is possible 
% to reconstruct the spectral reflectance bases s (and thus the full
% spectral reflectance S)

clear, clc, close all

%% For illuminant and reflectance, load basis functions and define:
% 1. Basis functions
% 2. Basis weights
% 3. Full spectra

plt_lin = 1; %want to plot?

load B_vrhel %basis functions of reflectance measurements
B_vrhel = B_vrhel(:,1:3); %Let's pretend it only has 3 basis functions
s = [-1;0.2;0.6]; %weights (arbitrary random picks)
S = B_vrhel*s; %Spectral reflectance function

load B_cieday %basis functions of daylight measurements
B_cieday = SplineSpd(S_cieday,B_cieday,S_vrhel); %interpolate to match sampling range and interval
S_cieday = S_vrhel;
[E,~,~,e(2),e(3)] = GenerateCIEDay(6500,B_cieday); %E = Spectral power distribution. e = weights.
e(1) = 1; %GenerateCIEDay sets e(1) = 1


if plt_lin %visualise the basis functions and the full spectra
    figure, hold on
    axis tight
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
end

%% Sensor Spectral Sensitivities

load T_cones_sp % Load sensor receptivities
T_cones_sp = SplineCmf(S_cones_sp,T_cones_sp,S_vrhel); %interpolate to match sampling range and interval
S_cones_sp = S_vrhel;
%plot(T_cones_sp')
R = T_cones_sp;

%%
% (eq numbers from second ref)

% Eq. 9.3
L = E.*S; %color signal = illum .* reflectance
p = R*L; %photon catch = sensor sensitivities * color signal

% % written alternatively:
p_alt1 = R*(E.*S); %photon catch = sensor sensitivities * illum .* reflectance

% Eq. 9.4
% Refers to intrisic colour calculation

% Eq. 9.5
if ~isequal(E, B_cieday*e')
    error('Error: SPD does not equal SPD reconstructed from basis functions')
end

% Eq. 9.6
if ~isequal(S, B_vrhel*s)
    error('Error: SRF does not equal SRF reconstructed from basis functions')
end

% Jumping back to thesis
% First part of Eq. at bottom of page 63:
p_alt2 = sum(E.*S.*R'); %with spectra
p_alt3 = sum(B_cieday*e'.*B_vrhel*s.*R'); %with linear functions and weights

% Second part

g=zeros(3,3,3); %preallocate
for i=1:3
    for j=1:3
        for k=1:3
            g(i,j,k) = sum(B_cieday(:,i).*B_vrhel(:,j).*R(k,:)');
        end
    end
end

%trace version
p_alt4 = [trace(s*e*g(:,:,1)),trace(s*e*g(:,:,2)),trace(s*e*g(:,:,3))];

% forloop version
p_alt5=[0;0;0];
for i=1:3
    for j=1:3
        for k=1:3
        p_alt5(k) = p_alt5(k) + s(j)*e(i)*g(i,j,k);
        end
    end
end

%would be nice to do this neatly in a single matrix multiplication

p_alt6(1) = e*g(:,:,1)*s;
p_alt6(2) = e*g(:,:,2)*s;
p_alt6(3) = e*g(:,:,3)*s;
%that's neat enough :)

%although:
eg = [e*g(:,:,1);e*g(:,:,2);e*g(:,:,3)]; %depends on light and stable matrix g
p_alt7 = eg*s;
%This seems equivalent to eq 9.7 (first ref) where the light and the
%independent matrix are bundled together, and it sure is neat.

%% Inverse

% Here I want to try to get the values for s (from which I can easily
% construct S) from p, e and g alone.

eg_inv = inv(eg); 
s_recovered = eg_inv*p;

% alternative formulation:
%s_recovered = eg\p;

s
s_recovered

%s-s_recovered % Essentially the same. Presumably just precision errors


