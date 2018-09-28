base = 'C:\Users\cege-user\Documents\Large data\Foster Images\';
for i=1:4 %2002 images
    ims(i)=load([base, '2002\scene',num2str(i),'.mat']); %imageS
end
%     %2004 images
ims(5)=load([base,'2004\scene1\ref_crown3bb_reg1.mat']);
ims(6)=load([base,'2004\scene2\ref_ruivaes1bb_reg1.mat']);
ims(7)=load([base,'2004\scene3\ref_mosteiro4bb_reg1.mat']);
ims(8)=load([base,'2004\scene4\ref_cyflower1bb_reg1.mat']);
ims(9)=load([base,'2004\scene5\ref_cbrufefields1bb_reg1.mat']);

%%

% This is a very rough fudge to get colour images. 
% 25,17,7 chosen without much thought.

% % for full images:
% for i = [5,6]% 1:length(ims)
%     clear disp
%     disp(:,:,1) = ims(i).reflectances(:,:,25);
%     disp(:,:,2) = ims(i).reflectances(:,:,17);
%     disp(:,:,3) = ims(i).reflectances(:,:,7);
%     
%     figure(i)
%     imshow(disp)
%     pause(1)
% end

% for cropping out the grey card:
for i = [5,6]% 1:length(ims)
    clear disp
    disp(:,:,1) = ims(i).reflectances(:,1:end-100,25); %roughly
    disp(:,:,2) = ims(i).reflectances(:,1:end-100,17);
    disp(:,:,3) = ims(i).reflectances(:,1:end-100,7);
    
    figure(i)
    imshow(disp)
    pause(1)
end

% I like 5 for mixed lighting and 6 for focal object in natural
% surroundings (and fairly even lighting I think

%%

S = [400,10,33];

%t = reshape(ims(5).reflectances,1018*1339,33);%temp
t = reshape(ims(6).reflectances,1017*1338,33);%temp
%t = reshape(ims(1).reflectances,748*820,31);%temp

% % visualise spectra
% figure, hold on
% for j=1:100
%     clf
%     plot(SToWls(S),t(j:50000:size(t,1),:))
%     drawnow
%     pause(0.1)
% end


%%

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(t);
figure, hold on
plot(COEFF(:,1:6))
legend('Location','best')

% Looks like removing the first and final few readings would drastically
% change the PCA, perhaps simply because high noise in recordings at these
% wavelengths
[COEFF2, SCORE2, LATENT2, TSQUARED2, EXPLAINED2] = pca(t(:,3:30));
figure, hold on
plot(COEFF2(:,1:6))
legend('Location','best')

% Would be nice to do this more neatly:
% Weight the data by a function that took cone sensitvities into account
% somehow. Or is that too circular?

% comparison of 1st PC to mean of data. Similar but not the same (I think:
% because pca.m uses svd as default)
figure, hold on
plot(COEFF(:,1)/max(COEFF(:,1)),'r--','DisplayName','normalised 1st PC')
plot(mean(t)/max(mean(t)),'g--','DisplayName','normalised mean')
legend('Location','best')


%% Correlation matrices

%t = reshape(ims(5).reflectances,1018*1339,33);%temp
t = reshape(ims(6).reflectances,1017*1338,33);%temp

c = corr(t);

figure
imagesc(c)
axis image
colorbar
colormap gray

S_refs_f = SToWls(S);
set(gca,'XTickLabel',S_refs_f(xticks))
set(gca,'YTickLabel',S_refs_f(xticks))

%% Diagramatic plot showing a comparison between s-cone sensitivity, 
% correlation between reflectances around peak s-cone sensitvity and other wavelengths, 
% and the second principal component of daylight

load B_cieday.mat
load T_cones_sp.mat

figure, hold on
plot(SToWls(S_cones_sp),T_cones_sp(3,:),'k--','LineWidth',3)
plot(SToWls(S), c(3:7,:))
plot(SToWls(S_cieday),B_cieday(:,2)/max(B_cieday(:,2)),'r:','LineWidth',3)




