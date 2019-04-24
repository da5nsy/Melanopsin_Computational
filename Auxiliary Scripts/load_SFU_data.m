function data = load_SFU_data(base)

% loads SFU data, preferentially from a .mat file, but if not found it will
% unzip the downloaded files, put all the data in a struct and save it out
% (for next time).

% Barnard, K., Martin, L., Funt, B. and Coath, A., 2002. A data set for color research. Color Research & Application, 27(3), pp.147–151.
% Barnard, K., Cardei, V. and Funt, B., 2002. A Comparison of Computational Color Constancy Algorithms - Part I: Methodology and Experiments With Synthesized Data. IEEE Transactions on Image Processing, 11(9), pp.972–984.

% http://www.cs.sfu.ca/~colour/data/colour_constancy_synthetic_test_data/index.html

% See also BarnardFunt_dataExtraction.m
% https://github.com/da5nsy/Melanopsin_Computational/blob/3b07be65789ff4228d36d5c4de474ad5cf52263f/Auxiliary%20Scripts/BarnardFunt_dataExtraction.m

%% Unzip data

if ~exist('base','var')
    base = 'C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\SFU\';
end

if exist([base,'data.mat'],'file') == 2
    load([base,'data.mat'],'data')
    return
end

gzFiles = dir([base,'*','gz']);

for i = 1:length(gzFiles)
gzFilesUnzipped(i) = gunzip([gzFiles(i).folder,'\',gzFiles(i).name]);
end

%%

for i = 1:length(gzFiles)
    f1 = fopen(gzFilesUnzipped{i});
    c = textscan(f1,'%s','Delimiter','\n');
    fclose(f1);
    data(i).data = str2double(c{1});
    data(i).data = data(i).data(~isnan(data(i).data));
    data(i).data = reshape(data(i).data,101,[]);
    %figure, plot(data)
    data(i).name = gzFiles(i).name(1:end-3);
end

save('data.mat','data')

end