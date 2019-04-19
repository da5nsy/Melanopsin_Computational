clear, clc, close all

[NUM,TXT,RAW] = xlsread('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Natural Reflectances\natural400_700_5.csv');

%%
clc

for surf = 1:218
    clear A B C
    range = 6*surf-4:6*surf-1;
    %range = 2:5;
    
    for i = 1:length(range)
        A = char(RAW{range(i),1});
        B = regexp(A,'\d*','Match');
        C(i,:) = str2double(B);
    end
    C = C';
    C = C(:);
    C(end+1) = RAW{6*surf,1};
    
    SRF(:,surf) = C;
    
end

figure, plot(SRF) %that does not look right, but at least I'm getting at the data

%6x-4 gives the first line



%% Things that should work but don't:

%data = importdata('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Natural Reflectances\natural400_700_5.asc')
%data = readmatrix('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Natural Reflectances\natural400_700_5.asc');










