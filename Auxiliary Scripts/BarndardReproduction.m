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

%nSurfList = [4, 8, 16, 32, 65, 128, 256, 512, 1024]; % Barnard et al.(10.1109/TIP.2002.802531)
nSurfList = [2, 4, 8, 16, 32, 64]; % Hordley & Finlayson (10.1364/JOSAA.23.001008)
nScenes = 1000;

%%

rng(1); %fixes random number generator for reproducibility

mink_norm = [-1,1]; %This switches the cc algo from max-RGB to grey-world
AE = NaN(length(mink_norm),length(nSurfList));
for mn = 1:length(mink_norm)
    for nSurfc = 1:length(nSurfList) %NSurfCount
        nSurf = nSurfList(nSurfc);
        AngularError = [];
        for scene = 1:nScenes
            
            % Generate random surfaces and illuminants and compute resulting RGBs
            % note: randi is sampling with replacement. Use randperm if you require all unique values.
            surfs = data(3).data(:,randi(size(data(3).data,2),nSurf,1));
            ill = data(5).data(:,randi(size(data(5).data,2)));
            RGB = data(4).data'*(ill.*surfs);
            ill_RGB = data(4).data'*ill;
            
            % RGB(1,:) = RGB(1,:)+1000; %artifically make correction more extreme - for testing only
            
            % Copmute algo on RGB
            [white_R ,white_G ,white_B] = general_cc(RGB,0,mink_norm(mn),0); %grey world
            
            %assess algo
            AngularError = [AngularError, atan2d(norm(cross(ill_RGB,[white_R ,white_G ,white_B])),dot(ill_RGB,[white_R ,white_G ,white_B]))];
            %rgDist = sqrt(
        end
        %calculate average performance for given number of surfs
        %average(number)
        AE(mn,nSurfc) = rms(AngularError); %Root meas squre Angular Error
    end
end

figure, hold on
plot(log2(nSurfList),AE(1,:),'o-','Color',[1 0.4 0.9],'MarkerFaceColor',[1 0.4 0.6],'DisplayName','Max-RGB')
plot(log2(nSurfList),AE(2,:),'gs-','MarkerFaceColor','g','DisplayName','Grey world')
xticks(log2(nSurfList))
xlabel('log 2 Num.Surfaces per Image')
ylabel('RMS angular error')
legend


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

