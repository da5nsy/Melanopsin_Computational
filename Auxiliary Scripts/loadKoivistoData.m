% Load and save Koivisto data

% Available from http://www.uef.fi/web/spectral/natural-colors (2019/04/25)
% Reference: Parkkinen, J., Jaaskelainen, T. and Kuittinen, M., 1988. Spectral representation of color images. In: [1988 Proceedings] 9th International Conference on Pattern Recognition. [online] [1988 Proceedings] 9th International Conference on Pattern Recognition. Rome, Italy: IEEE Comput. Soc. Press, pp.933–935. Available at: <http://ieeexplore.ieee.org/document/28405/> [Accessed 5 Oct. 2018].

%%
clear, clc, close all

% Load data (converted to csv previously)
[NUM,TXT,RAW] = xlsread('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Natural Reflectances\natural400_700_5.csv');

%%

for surf = 1:219
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
    
    sur_koivisto(:,surf) = C;
    labels_koivisto(surf).label = RAW{range(i)-4,1};
    
end

sur_koivisto(sur_koivisto > 4096) = 4096; % As per note on website

S_koivisto = WlsToS([400:5:700]');

figure, plot(SToWls(S_koivisto),sur_koivisto)

%% Save in PTB format for further use

save('sur_koivisto.mat','sur_koivisto','labels_koivisto','S_koivisto')

%% Check it all worked:

% clear, clc, close all
% 
% load sur_koivisto.mat
% 
% figure, plot(SToWls(S_koivisto),sur_koivisto)











