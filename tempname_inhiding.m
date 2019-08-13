function [output,mark,sf_l,sf_s] = tempname_inhiding(lsri,Lum,autoScale,pcSurfRange)

%%

autoScale = 0;
pcSurfRange = 20:10:100;

%%
clc

nMethods = 4;

nSurf       = zeros(length(pcSurfRange),1);
sel_store   = zeros(length(pcSurfRange),round((pcSurfRange(end)/100) * size(lsri,2)),size(lsri,3));
output      = zeros(2,size(lsri,2),size(lsri,3),nMethods,length(pcSurfRange));
sf_l        = zeros(length(pcSurfRange),1);
sf_s        = zeros(length(pcSurfRange),1);
mark        = zeros(nMethods,length(pcSurfRange));

for pcSI = 1:length(pcSurfRange) %percent surfaces index
    nSurf(pcSI) = round((pcSurfRange(pcSI)/100) * size(lsri,2)); %number of surfaces
    sel   = zeros(nSurf(pcSI),size(lsri,3));
    lsri2 = zeros(size(lsri,1),nSurf(pcSI),size(lsri,3));
    Lum2  = zeros(nSurf(pcSI),size(lsri,3));
    for i = 1:size(lsri,3)
        sel(:,i)     = randperm(size(lsri,2),nSurf(pcSI)); %selection
        lsri2(:,:,i) = lsri(:,sel(:,i),i);
        Lum2(:,i)    = Lum(sel(:,i),i);
    end
    sel_store(pcSI,1:nSurf(pcSI),:) = sel;
    
    % Perform corrections
    [output(:,1:nSurf(pcSI),:,:,pcSI),sf_l(pcSI),sf_s(pcSI)] = performCC(lsri2,Lum2,autoScale);
   
    % Score corrections
    for i=1:nMethods
        rng(1)
        mark(i,pcSI) = KMeansMark(squeeze(output(:,1:nSurf(pcSI),:,i,pcSI)),size(T_SRF,2),sel);
    end
    
    disp(pcSurfRange(pcSI))
    
    %     plot(pcSurfRange,mark','.')
    %     xlim([0 100])
    %     ylim([0 1])
    %     drawnow
end
end