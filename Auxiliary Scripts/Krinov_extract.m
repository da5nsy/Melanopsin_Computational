%% Extracting Krinov data from data provided by David Brainard.

% The Canadian translation of Krinov's work (http://doi.org/10.4224/20386770) 
% lists 370 measurements (though a small number are arguably non-natural).

% Data provided by David Brainard provides 191 spectra (2 of which are differently scaled).

% Description from Kohonen+ doi.org/10.1002/col.20244 :
% "The Krinov data consists of 337 surface reflectance spectra
% acquired from the naturally occurring objects such as forests,
% shrubs, roads, grass, snow, and water surfaces. The
% spectral range of the original data was from 400 to 650
% nm at 10 nm intervals. The used data was downloaded
% from http://vision.cs.anzona.edu/kobus/research/data/colour_
% constancy_synthetic_test_data/, and it was converted to
% match the output format of the PhotoResearch PR-630 spectrophotometer,
% from 380 to 780 nm at 4 nm intervals. The
% Krinov data has been used especially for early investigations
% of the effective dimension of spectral reflectances."

% The SFU dataset mentioned above
% (now at http://www.cs.sfu.ca/~colour/data/colour_constancy_synthetic_test_data/)
% however, says it contains 350 Krinov surfaces (though appears to actually
% include 355...)
% Extraction script here: https://github.com/da5nsy/Melanopsin_Computational/blob/2d324882e7461613c197fdd52b38951e9865224a/Auxiliary%20Scripts/BarnardFunt_dataExtraction.m


%%

clear, clc, close all

dataFileName = 'C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Brainard_Krinov\krinov.text';
fid = fopen(dataFileName, 'r');
sur_krinov = fscanf(fid,'%f',[191,Inf])';
fclose(fid);

sur_krinov(:,14) = sur_krinov(:,14)/1000;
sur_krinov(:,186) = sur_krinov(:,186)/1000;

S_krinov = WlsToS([400:10:650]');

% figure,
% plot(SToWls(S_krinov),sur_krinov');

%%

clear key labels_krinov
clc

keyFileName = 'C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Brainard_Krinov\krinov.key';

fid = fopen(keyFileName, 'r');
while true
    thisline = fgetl(fid);
    if exist('key')
        key{end+1} = thisline;
    else
        key = {thisline};
    end
    if ~ischar(thisline); break; end  
end
fclose(fid);
key = key(:,1:end-1)';

labels_krinov = cell2struct(key','label')';

%%

clear thisline fid ans dataFileName key keyFileName

save('sur_krinov.mat')

