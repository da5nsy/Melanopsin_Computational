function [output,sf_l,sf_s] = performCC(lsri,Lum,autoScale)

% Do nothing 'correction' 
% -----------------------
% Baseline measurement
output(:,:,:,1) = lsri(1:2,:,:);


% Grey-world assumption:
% ----------------------
% Here, the average chromaticity of the reflectances under a single 
% illuminant will be taken to be the origin, and everything will be denoted
% in reference to that.

gw = zeros(size(lsri(1:2,:,:)));
for i=1:size(lsri,3)
    IllEst = mean(lsri(1:2,:,i),2); %Illuminant estimate
    gw(:,:,i) = lsri(1:2,:,i) - IllEst;
end
output(:,:,:,2) = gw;

% % Figure to check operation
% figure, hold on
% scatter(lsri(1,:,end),lsri(2,:,end),'k')
% scatter(IllEst(1),IllEst(2),'r','filled')
% scatter(gw(1,:,end),gw(2,:,end),'ks')
% plot([0,0],[min(ylim),max(ylim)],'k--')
% plot([min(xlim),max(xlim)],[0,0],'k--')
% scatter(0,0,'rs','filled')


% Bright-is-white assumption:
% ---------------------------
% Here, the surface with the highest luminance value will be taken as being
% a neutral chromaticity, and will be shifted to the origin, with
% everything else denoted in reference to that.

BiW = zeros(size(lsri(1:2,:,:)));
maxLumLoc = zeros(size(lsri,3),1);
for i=1:size(lsri,3)
    [~,maxLumLoc(i)] = max(Lum(:,i));
    IllEst(:,i) = lsri(1:2,maxLumLoc(i),i); %Illuminant estimate
    BiW(:,:,i) = lsri(1:2,:,i) - IllEst(:,i);
end
output(:,:,:,3) = BiW;

% % Figures to check operation
% figure, hold on
% for i=1:500:size(lsri,3)
% scatter(lsri(1,:,i),lsri(2,:,i))
% ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex-1; %resets colour order so that colours correspond
% scatter(BiW(1,:,i),BiW(2,:,i),'s')
% end
% plot([0,0],[min(ylim),max(ylim)],'k--')
% plot([min(xlim),max(xlim)],[0,0],'k--')
% figure, 
% hist(maxLumLoc,100)


% Melanopsin based correction
% ---------------------------

if autoScale == 1
    [sf_l,sf_s] = melcomp_6_calcsf(lsri, -2:0.01:2,-2:0.01:2); %calculates scaling factors
else
    sf_l = 1.081;
    sf_s = -0.797;
end

MC = lsri;
MC(1,:) = MC(1,:)+sf_l*MC(4,:);
MC(2,:) = MC(2,:)+sf_s*MC(4,:);
output(:,:,:,4) = MC(1:2,:,:);

end