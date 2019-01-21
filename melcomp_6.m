clear
clc

%% Load data

% Tables from CIE 170-2:2015

% live:     http://files.cie.co.at/813_Tables_CIE_170-2.xls
% archive:  https://web.archive.org/web/20190121134312/http://files.cie.co.at/813_Tables_CIE_170-2.xls

T_CIE_MB_10 = xlsread('C:\Users\cege-user\Dropbox\UCL\Data\Colour standards\813_Tables_CIE_170-2.xls','Table 10.6','B2:D90');
S_CIE_MB_10 = [390,5,89];

%% Plot spectral locus

figure,
scatter(T_CIE_MB_10(:,1),T_CIE_MB_10(:,3))
xlim([0.4,1])

%% 

EE = ones(length(T_CIE_MB_10),1);

T_CIE_MB_10'*EE;


