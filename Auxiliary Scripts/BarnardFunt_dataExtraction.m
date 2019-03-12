% Code to read Barnard/Funt Reflectance data

% http://www2.cs.sfu.ca/~colour/data/colour_constancy_synthetic_test_data/index.html
% http://www.cs.sfu.ca/~colour/data/colour_constancy_synthetic_test_data/reflect_db.reflect.gz

% "The surface reflectance data (reflect_db.reflect) is a set of 1995 spectra 
% compiled from several sources. These surfaces include 
% the 24 Macbeth color checker patches, 
% 1269 Munsell chips, 
% 120 Dupont paint chips [1], 
% 170 natural objects [1], 
% the 350 surfaces in Krinov data set [2], and 
% 57 additional surfaces measured by ourselves."

%24+1269+120+170+350+57 =1990 (not 1995...)

clc, clear, close all

%first, rename as .csv, open in excel and save as an xlsx doc.
%I'm sure there's a better way to do this but this works (slowly)

%M = xlsread('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\reflect_db.reflect.xlsx');

%save...?

%%

load('reflect_db.mat') %if you saved previously

idx = all(isnan(M),2);
idr = diff(find([1;diff(idx);1]));
D = mat2cell(M,idr(:),size(M,2));
D = D(1:2:end);

%%

D = cell2mat(D');

%working out what/where everything is
%new styles seem to start at:
%(1) 24+57 (Macbeth +own)
%82 %120 (DuPont)
%202 %350(355) (Krinov) %!!!!!!!!!!!!!!!!
%557 %1269 (Munsell)
%1826
%(end 1995)

%%

%I've got the data just about in shape, but I'd be cautious about using
%them since I don't know particularly well what the sampling range/interval
%is or what each material is. (Or why there are 5 more Krinov measurements
%than there should be. Though - there are 370 in the Krinov publication, so
%who knows...)
