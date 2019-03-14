function [FMS] = ForsythMeasurementOfSuccess(lsri,pltc_alt)

% Forsyth, D.A., 1990. A novel algorithm for color constancy. 
% International Journal of Computer Vision, 5(1), pp.5–35.
% https://doi.org/10.1007/BF00056770


%%  median distance of the outputs from the average 

%average
av = mean(lsri,3); %average chromaticity per object
figure, hold on
scatter(lsri(1,:),lsri(2,:),...
    [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
scatter(av(1,:),av(2,:),'k.')

%distances from average
D = zeros([size(lsri,2),size(lsri,3)]);
for i=1:size(lsri,2)
    for j= 1:size(lsri,3)
        plot([av(1,i),lsri(1,i,j)],[av(2,i),lsri(2,i,j)],'k')
        D(i,j) = sqrt((av(1,i)-lsri(1,i,j))^2 + (av(2,i)-lsri(2,i,j))^2); %distance
    end
end

Dm = median(D,2); %distance median

%% normalised by the euclidian magnitude of the outputs

em = sqrt(sum(D.^2,2)); %euclidian magnitude

FMS = mean(Dm./em);


