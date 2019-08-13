function [output,mark,sf_l,sf_s,sel_store] = performCCandMarkForRange(lsri,Lum,autoScale,pcSurfRange)

nMethods = 4;

nSurf       = zeros(length(pcSurfRange),1);
sel_store   = zeros(length(pcSurfRange),round((pcSurfRange(end)/100) * size(lsri,2)),size(lsri,3));
output      = zeros(2,length(pcSurfRange),size(lsri,3),nMethods,length(pcSurfRange));
sf_l        = zeros(length(pcSurfRange),1);
sf_s        = zeros(length(pcSurfRange),1);
mark        = zeros(nMethods,length(pcSurfRange));

for pcSI = 1:length(pcSurfRange) %percent surfaces index
    nSurf(pcSI) = round((pcSurfRange(pcSI)/100) * size(lsri,2)); %number of surfaces
    sel   = zeros(nSurf(pcSI),size(lsri,3));
    lsri2 = zeros(size(lsri,1),nSurf(pcSI),size(lsri,3));
    Lum2  = zeros(nSurf(pcSI),size(lsri,3));
    for ill = 1:size(lsri,3)
        sel(:,ill)     = randperm(size(lsri,2),nSurf(pcSI)); %selection
        lsri2(:,:,ill) = lsri(:,sel(:,ill),ill);
        Lum2(:,ill)    = Lum(sel(:,ill),ill);
    end
    sel_store(pcSI,1:nSurf(pcSI),:) = sel;
    
    % Perform corrections
    [output(:,1:nSurf(pcSI),:,:,pcSI),sf_l(pcSI),sf_s(pcSI)] = performCC(lsri2,Lum2,autoScale,sel);
    
    % Score corrections
    for meth=1:nMethods
        %rng(1)
        mark(meth,pcSI) = KMeansMark(squeeze(output(:,1:nSurf(pcSI),:,meth,pcSI)),size(lsri,2),sel);
    end
    
    disp(pcSurfRange(pcSI))
    
    %     plot(pcSurfRange,mark','.')
    %     xlim([0 100])
    %     ylim([0 1])
    %     drawnow
end
end