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

for i = [5,6]% 1:length(ims)
    clear disp
    disp(:,:,1) = ims(i).reflectances(:,:,25);
    disp(:,:,2) = ims(i).reflectances(:,:,17);
    disp(:,:,3) = ims(i).reflectances(:,:,7);
    
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

figure, hold on
for j=1:100
    clf
    plot(SToWls(S),t(j:50000:size(t,1),:))
    drawnow
    pause(0.1)
end


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




