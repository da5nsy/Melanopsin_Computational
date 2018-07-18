clear, clc, close all

% A fresh attempt

%% Load Image

load sur_vrhel
refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
sur_vrhel=sur_vrhel(:,refs);

% base = 'C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Foster Lab Images\';
% im=load([base, '2002\scene',num2str(1),'.mat']); 
% %im = im.reflectances(1:25,1:25,:);
% im = im.reflectances(1:748,:,:);
% S_im=[410,10,31];
% %figure,imagesc(im.reflectances(:,:,17)); colormap('gray'); brighten(0.5);

%% Load Daylight SPD
D_CCT=1./linspace(1/3600,1/25000,20); %non-linear range, aiming to better reproduce observed variation
load B_cieday
daylight_spd = GenerateCIEDay(D_CCT,[B_cieday]); %these appear to be linearly upsampled from 10nm intervals
daylight_spd = daylight_spd./repmat(max(daylight_spd),81,1); %normalise
%figure, plot(daylight_spd)

%% Observer

% Smith-Pokorny, for use with MacLeod Boynton diagram
load T_cones_sp
load T_rods
load T_melanopsin

% Pull them all together
T_LMSRI=[(SplineCmf(S_cones_ss2,T_cones_ss2,S_im));...
    (SplineCmf(S_rods,T_rods,S_im));...
    (SplineCmf(S_melanopsin,T_melanopsin,S_im))];
S_LMSRI=S_im;

%figure, plot(SToWls(S_LMSRI),T_LMSRI)

%% Convert to radiance
spd_i=SplineSpd(S_cieday,daylight_spd,S_im,1); %interpolate to match range and interval of Foster images

im_r=zeros([size(im),length(D_CCT)]);

for i = 1:size(im,3)
    for j=1:length(D_CCT)
        im_r(:,:,i,j) = im(:,:,i)*spd_i(i,j); %image radiance
    end
end

%% Calculate LMSRI and lsri for each pixel
[r,c,w,d] = size(im_r);

im_LMSRI = zeros(r*c,5,d); %First level
im_lsri  = zeros(r*c,4,d); %Second level
im_rr    = zeros(r*c,31,d);

for i=1:20
    im_rr(:,:,i) = reshape(im_r(:,:,:,i), r*c, w);

    im_LMSRI(:,:,i) = (T_LMSRI*im_rr(:,:,i)')'; 
    im_lsri(:,:,i)  = im_LMSRI(:,[1,3,4,5],i)./(im_LMSRI(:,1,i)+im_LMSRI(:,2,i));
end

%%
figure, hold on
for i=1:20
    scatter3(im_lsri(:,1,i),im_lsri(:,2,i),ones(size(im_lsri,1), 1)*i,'filled')
    drawnow
end
%axis equal

%%

im_LMSRI_c = log(im_LMSRI)-mean(log(im_LMSRI)); %'c' = 'corrected'
im_lsri_c  = log(im_lsri)-mean(log(im_lsri));   %'c' = 'corrected'

figure, hold on
for i=1:20
    scatter3(im_lsri_c(:,1,i),im_lsri_c(:,2,i),ones(size(im_lsri,1), 1)*i,'filled')
    drawnow
end

%%

figure, hold on
for i=1:5:625
    plot3(squeeze(im_lsri(i,1,:)),squeeze(im_lsri(i,2,:)),1:20)
end

% figure, hold on
% for i=1:5:625
%     plot3(squeeze(im_lsri_c(i,1,:)),squeeze(im_lsri_c(i,2,:)),1:20)
% end

%%
figure, hold on
for i=1:25
    for j=1:25
        plot(SToWls(S_im),squeeze(im(i,j,:))/max(squeeze(im(i,j,:))))
    end
end

%%

im_r2=reshape(im,625,31);

% figure, hold on
% for i=1:625
%     plot(SToWls(S_im),im_r2(i,:))
% end

figure, hold on
s=std(im_r2);
plot(SToWls(S_im),s)

%%

%%

im_r2=reshape(im,748*820,31);

% figure, hold on
% for i=1:625
%     plot(SToWls(S_im),im_r2(i,:))
% end

figure, hold on
s=std(im_r2);
plot(SToWls(S_im),s)

%%
clear, clc, close all

base = 'C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Foster Lab Images\';
for i=1:4 %2002 images
    ims(i)=load([base, '2002\scene',num2str(i),'.mat']); %imageS
end
%2004 images
ims(5)=load([base,'2004\scene1\scene1\ref_crown3bb_reg1.mat']);
ims(6)=load([base,'2004\scene2\scene2\ref_ruivaes1bb_reg1.mat']);
ims(7)=load([base,'2004\scene3\scene3\ref_mosteiro4bb_reg1.mat']);
ims(8)=load([base,'2004\scene4\scene4\ref_cyflower1bb_reg1.mat']);
ims(9)=load([base,'2004\scene5\scene5\ref_cbrufefields1bb_reg1.mat']);

ims(1).reflectances = ims(1).reflectances(1:748,:,:);
ims(2).reflectances = ims(2).reflectances(1:700,:,:);
ims(3).reflectances = ims(3).reflectances(1:750,:,:);
ims(4).reflectances = ims(4).reflectances(1:664,96:end,:);

figure, hold on
for i=1:9
    ims(i).reshape=log(reshape(ims(i).reflectances,size(ims(i).reflectances,1)*size(ims(i).reflectances,2),size(ims(i).reflectances,3)));
    
    if size(ims(i).reflectances,3)==31 %2002 images
        ims(i).S = [410,10,31];
    elseif size(ims(i).reflectances,3)==33 %2004 images #1:4
        ims(i).S = [400,10,33];
    elseif size(ims(i).reflectances,3)==32 %2004 image #5
        ims(i).S = [400,10,32];
    end
    
    ims(i).std=nanstd(ims(i).reshape);
    plot(SToWls(ims(i).S),ims(i).std/max(ims(i).std))
end

figure, hold on
for i=1:9
    ims(i).m=nanmean(ims(i).reshape); %m for mean
    plot(SToWls(ims(i).S),ims(i).std./ims(i).m)
end
grid on


