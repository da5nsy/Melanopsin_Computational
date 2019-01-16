function Mel_Phot_Corr()

% Is there a greater correlation between melanopic and photopic luminance
% for natural light sources than for artificial light sources?

%% Pre-flight

clear, clc, close all

PF_justReal = 1; %Exclude theoretical illums from Houser data

%% Load Data

% Load artifical lighting spds
% https://github.com/Psychtoolbox-3/Psychtoolbox-3/blob/0dac597065c348d0efb8d024c052fd8e9ff17322/Psychtoolbox/PsychColorimetricData/PsychColorimetricMatFiles/spd_houser.mat

% 401 normalised illuminant spectral power distributions from:
% "Review of measures for light-source color rendition and considerations
% for a two-measure system for characterizing color rendition"
% Kevin W. Houser, Minchen Wei, Aurélien David, Michael R. Krames, and Xiangyou Sharon Shen
% Optics Express, Vol. 21, Issue 8, pp. 10393-10411 (2013)
% http://dx.doi.org/10.1364/OE.21.010393

load spd_houser.mat

codes_houser = strsplit('H	H	G	C	C	D	E	E	E	E	F	F	H	H	H	H	H	H	H	H	H	A	A	G	G	G	G	L	L	L	L	C	C	C	C	C	C	C	C	C	D	D	D	E	E	D	D	D	D	D	C	C	C	C	C	C	C	C	C	C	C	E	E	E	E	E	H	H	H	H	H	H	I	I	I	I	H	H	H	H	B	B	B	B	B	H	H	H	H	G	H	H	H	H	H	H	H	H	H	G	G	G	G	G	G	G	G	G	G	G	G	G	G	G	G	G	A	A	A	A	A	A	A	J	J	J	J	J	J	J	J	K	K	K	K	K	K	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	E	H	H	H	H	H	H	G	G	G	G	G	H	H	H	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	I	E	E	E	E	F	H	B	A	A	B	E	E	E	F	D	B	B	B	B	B	B	B	I	H	H	H	H	H	H	H	F	F	F	D	D	B	B	B	D	D	F	D	D	D	D	D	D	D	D	D	D	D	D	C	C	F	F	F	F	E	E	E	E	L	L	D	D	D	D	D	C	C	C	C	C	C	F	E	E	E	F	E	E	E	E	A	A	A	A	A	A	A	A	A	A	A	A	A	A	A	E	F	F	F	A	A	A	A	G');
% could do above with xlsread, would be neater but slower, and relies on
% others using this script having the xls
 
% 1-Letter Code	Lamp Category	
% 	Real Illuminants	
% A		LED Phosphor Real
% B		LED Mixed Real
% C		Fluorescent Broadband
% D		Fluorescent Narrowband
% E		HID
% F		Tungesten Filament
% 	Theoretical Illuminants	
% G		LED Phosphor Models
% H		LED Mixed Models
% I		Fluorescent Models
% J		Blackbody Radiation
% K		D-Series Illuminant
% L		Other (e.g., Equal-Energy, Clipped Incan, Ideal Prime Color)

if PF_justReal
    spd_houser = spd_houser(:,ismember(codes_houser,{'A','B','C','D','E','F'}));
end

T_artif = spd_houser;



% Load Daylight Data
load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
granada = final; clear final
% 300 - 1100nm, 5nm interval, unlabeled
% 2600 samples
T_day = granada(17:97,:); %match obs
S_day = [380,5,81];
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)

%% Load observer

% Obs data
load('T_cones_ss10')  % 10 deg obs
load('T_melanopsin')

%% Calculate cone/mel values

% Calculate LMS of daylight samples
LMS_d(:,:)=SplineCmf(S_cones_ss10,T_cones_ss10,S_day)*T_day;
LMS_a(:,:)=SplineCmf(S_cones_ss10,T_cones_ss10,S_day)*T_artif;

% Calulate M of daylight samples
Mel_d(:,:)=SplineCmf(S_melanopsin,T_melanopsin,S_day)*T_day;
Mel_a(:,:)=SplineCmf(S_melanopsin,T_melanopsin,S_day)*T_artif;

R_d = (LMS_d(1,:)+LMS_d(2,:))./Mel_d; %r for ratio
R_a = (LMS_a(1,:)+LMS_a(2,:))./Mel_a; 


%% 

% %lin x-axis
% figure, hold on
% histogram(R_d,'BinWidth',0.1,'Normalization','probability')
% histogram(R_a,'BinWidth',0.1,'Normalization','probability')
% xlabel('(L+M) / I')
% ylabel('Probability')

%log x-axis
R_dl=log10(R_d);
R_al=log10(R_a);
figure, hold on
histogram(R_dl,'BinWidth',0.02,'Normalization','probability')
histogram(R_al,'BinWidth',0.02,'Normalization','probability')
xlabel('log((L+M) / I)')
ylabel('Probability')

legend({'Daylight','Artificial'},'Location','best')


end