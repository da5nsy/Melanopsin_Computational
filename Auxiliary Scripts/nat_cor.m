% Correlation matrices of natural reflectances

clear, clc, close all

plot_where = [500,200];
plot_size  = [800,400];

base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\nat_corr';

p = 1; % 0 = figures display but don't save. 1 = figures display and save.
set(0,'defaultAxesFontName', 'Courier')

%% Load data

base = 'C:\Users\cege-user\Documents\Large data\Foster Images\';
for i=1:4 %2002 images
    ims(i)=load([base, '2002\scene',num2str(i),'.mat']); %imageS
end
%2004 images
ims(5)=load([base,'2004\scene1\ref_crown3bb_reg1.mat']);
ims(6)=load([base,'2004\scene2\ref_ruivaes1bb_reg1.mat']);
ims(7)=load([base,'2004\scene3\ref_mosteiro4bb_reg1.mat']);
ims(8)=load([base,'2004\scene4\ref_cyflower1bb_reg1.mat']);
ims(9)=load([base,'2004\scene5\ref_cbrufefields1bb_reg1.mat']);
%rename the field for clarity later
[ims.r] = ims.reflectances;
ims = rmfield(ims,'reflectances');

for i=1:length(ims) %no logic for "t". It was temp, and now it is permenent
    
    %reshape
    ims(i).t = reshape(ims(i).r,size(ims(i).r,1)*size(ims(i).r,2),size(ims(i).r,3));
    
    %add wavelength sampling info in PTB format
    if ismember(i,[1,2,3,4])
        ims(i).S = [410,10,31];
    elseif ismember(i,[5,6,7,8])
        ims(i).S = [400,10,33];
    elseif i==9
        ims(i).S = [400,10,32];
    end
    
end

load sur_vrhel.mat

% Not natural (?), but useful for comparison
% Not sure why the deadspace is included. Presumably so that they match but
% that seems like a silly idea.
load sur_nickerson.mat
sur_nickerson = sur_nickerson(5:65,:);
S_nickerson = [400,5,61];
load sur_macbeth.mat
sur_macbeth = sur_macbeth(2:71,:);
S_macbeth = [385,5,70];

%%
%close all

%visualization of corrs
figure; hold on
for i=1:length(ims)
    %calculate corr for each dataset
    ims(i).c = corr(ims(i).t);
    
    subplot(3,3,i)
    imagesc(ims(i).c)
    axis image
    colorbar
    colormap gray
    
    S_refs_f = SToWls(ims(i).S);
    set(gca,'XTickLabel',S_refs_f(xticks))
    set(gca,'YTickLabel',S_refs_f(xticks))
    
end

%visualization of images
figure; hold on
for i=1:length(ims)
    subplot(3,3,i)
    imshow(cat(3,ims(i).r(:,:,25),ims(i).r(:,:,17),ims(i).r(:,:,7)).^0.2) %added the final term as a bodge to make everything visible rather than for a strict scientific reason
end

figure, 
%median vs mean? Arguments for either, results similar but not same
av = median(cat(3,...
    ims(1).c,...
    ims(2).c,...
    ims(3).c,...
    ims(4).c,...
    ims(5).c(2:end-1,2:end-1),...
    ims(6).c(2:end-1,2:end-1),...
    ims(7).c(2:end-1,2:end-1),...
    ims(8).c(2:end-1,2:end-1),...
    ims(9).c(2:end,2:end)),...
    3);
imagesc(av)
axis image
colorbar
colormap gray

S_refs_f = SToWls(ims(1).S);
set(gca,'XTickLabel',S_refs_f(xticks))
set(gca,'YTickLabel',S_refs_f(xticks))

%% Plot all other datasets

%Something is not right with the axes labelling currently !!!!!!!!!

figure,

subplot(2,2,1)
sur_vrhel_c = corr(sur_vrhel');
imagesc(sur_vrhel_c)
axis image
colormap gray
S_refs_f = SToWls(S_vrhel);
set(gca,'XTickLabel',S_refs_f(xticks))
set(gca,'YTickLabel',S_refs_f(xticks))
xlabel('Vrhel')

subplot(2,2,2)
%sur_vrhel_n = sur_vrhel(:,[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141]);
sur_vrhel_n = sur_vrhel(:,[1:44,65,69,81:154]);
sur_vrhel_n_c = corr(sur_vrhel_n');
imagesc(sur_vrhel_n_c)
axis image
colormap gray
S_refs_f = SToWls(S_vrhel);
set(gca,'XTickLabel',S_refs_f(xticks))
set(gca,'YTickLabel',S_refs_f(xticks))
xlabel('Vrhel (Nat only)')

% There's a note in the PTB contents.m file that says that these are
% innacurate and in need of updating.
subplot(2,2,3)
sur_macbeth_c = corr(sur_macbeth');
imagesc(sur_macbeth_c)
axis image
colormap gray
S_refs_f = SToWls(S_macbeth);
set(gca,'XTickLabel',S_refs_f(xticks))
set(gca,'YTickLabel',S_refs_f(xticks))
xlabel('Macbeth')

subplot(2,2,4)
sur_nickerson_c = corr(sur_nickerson');
imagesc(sur_nickerson_c)
axis image
colormap gray
S_refs_f = SToWls(S_nickerson);
set(gca,'XTickLabel',S_refs_f(xticks))
set(gca,'YTickLabel',S_refs_f(xticks))
xlabel('Nickerson')

%%
% load T_cones_sp
% load T_melanopsin
% load T_rods
% figure, hold on
% subplot(2,2,1)
% plot(SToWls(S_cones_sp),T_cones_sp)
% plot(SToWls(S_melanopsin),T_melanopsin)
% plot(SToWls(S_rods),T_rods)
    
%% What about daylight SPDs?

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)
T_SPD=final; clear final
S_SPD=[300,5,161];

T_SPD(T_SPD < 0) = 0; 

T_SPD_c = corr(T_SPD');
T_SPD_c_log = corr(log10(T_SPD'));

figure,
imagesc(T_SPD_c)
colorbar

figure,
imagesc(T_SPD_c_log)
colorbar

