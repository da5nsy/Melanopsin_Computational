clear, clc, close all

% A fresh attempt

%% Load Image
base = 'C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Foster Lab Images\';
im=load([base, '2002\scene',num2str(1),'.mat']); 
im.reflectances = im.reflectances(1:748,:,:);
%figure,imagesc(im.reflectances(:,:,17)); colormap('gray'); brighten(0.5);

%% Load Daylight SPD
D_CCT=1./linspace(1/3600,1/25000,20); 
load B_cieday
daylight_spd = GenerateCIEDay(D_CCT,[B_cieday]); %these appear to be linearly upsampled from 10nm intervals
daylight_spd = daylight_spd./repmat(max(daylight_spd),81,1); %normalise
%figure, plot(daylight_spd)

%% Convert to radiance
spd_i=SplineSpd(S_cieday,daylight_spd,S_im,1); %interpolate to match range and interval of Foster images

% figure, hold on, %check interpolation
% scatter(SToWls(S_cieday),  spd)
% scatter(SToWls(S_im),      spd_i);

for i = 1:size(im,3)
    im_r(:,:,i) = im(:,:,i)*spd_i(i); %image radiance
end
