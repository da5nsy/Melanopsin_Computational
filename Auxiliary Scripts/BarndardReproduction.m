% Reproduce some of the results of:
% Barnard, K., Cardei, V. and Funt, B., 2002. A Comparison of Computational Color Constancy Algorithms - Part I: Methodology and Experiments With Synthesized Data. IEEE Transactions on Image Processing, 11(9), pp.972–984.

clear, clc, close all

%% Load data

data = load_SFU_data;
S_data = [380,4,101];

% for i = 1:length(data)
% figure,
% plot(SToWls(S_data),data(i).data)
% title(data(i).name, 'Interpreter', 'none')
% end

%%

nSurfList = [4, 8, 16, 32, 65, 128, 256, 512, 1024];
nScenes = 1000;

%%

rng(1); %fixes random number generator for reproducibility

for nSurfc = 1:length(nSurfList) %NSurfCount
    nSurf = nSurfList(nSurfc);
    for scene = 1:nScenes 
        
        % Generate random surfaces and illuminants and compute resulting RGBs
        % note: randi is sampling with replacement. Use randperm if you require all unique values.
        surfs = data(3).data(:,randi(size(data(3).data,2),nSurf,1));        
        ill = data(5).data(:,randi(size(data(5).data,2)));
        RGB = data(4).data'*(ill.*surfs);
        
        % RGB(1,:) = RGB(1,:)+1000; %artifically make correction more extreme - for testing only
        
        % Copmute algo on RGB
        %general_cc(RGB,0,1,0) %grey world
        %fixedRGB = x(RGB)
        
        %assess algo
        %number = metric(makeitRGB(ill),estimateWhitePoint(fixedRGB))
        %store number somewhere somehow
    end
    %calculate average performance for given number of surfs
    %average(number)
end




%%

[white_R ,white_G ,white_B,output_data] = general_cc(RGB,0,1,0);

figure, hold on
scatter3(RGB(1,:),RGB(2,:),RGB(3,:))
scatter3(output_data(1,:),output_data(2,:),output_data(3,:),'*')

%% Visualise correction

figure, 
subplot(1,2,1)
image1 = reshape(RGB,[3,32,32]);
image1 = permute(image1,[3,2,1]);
image1 = uint8(image1/max(image1(:))*255);
imshow(image1)

subplot(1,2,2)
image2 = reshape(output_data,[3,32,32]);
image2 = permute(image2,[3,2,1]);
image2 = uint8(image2/max(image2(:))*255);
imshow(image2)

