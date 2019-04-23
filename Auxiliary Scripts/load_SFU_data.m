%load_SFU_data

% http://www.cs.sfu.ca/~colour/data/colour_constancy_synthetic_test_data/index.html

% See also BarnardFunt_dataExtraction.m

%% Unzip data

base = 'C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\SFU\';
gzFiles = dir([base,'*','gz']);

for i = 1:length(gzFiles)
gzFilesUnzipped(i) = gunzip([gzFiles(i).folder,'\',gzFiles(i).name]);
end

%%

f1 = fopen(gzFilesUnzipped{1}); % your data in file test2.txt
c = textscan(f1,'%s','Delimiter','\n');
fclose(f1);
ImageDataSources_illum = str2double(c{1}(find(strcmp(c{1},'#'))+1:end));
ImageDataSources_illum = ImageDataSources_illum(~isnan(ImageDataSources_illum));
ImageDataSources_illum = reshape(ImageDataSources_illum,[101,11]);
figure, plot(ImageDataSources_illum)

%%
f1 = fopen(gzFilesUnzipped{2}); % your data in file test2.txt
c = textscan(f1,'%s','Delimiter','\n');
fclose(f1);
MeasuredWithSources_illum = str2double(c{1}(find(strcmp(c{1},'#! t=i n=101 o=380.0 s=4.0'))+1:end));
MeasuredWithSources_illum = MeasuredWithSources_illum(~isnan(MeasuredWithSources_illum));
MeasuredWithSources_illum = reshape(MeasuredWithSources_illum,[101,81]);
figure, plot(MeasuredWithSources_illum)

