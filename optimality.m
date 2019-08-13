function KMM = optimality(mel_offset)

%% Load Data

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_1',...
    'SSF','SS10',...
    'lum','CIE_10',...
    'mel_offset',mel_offset);

%% Colorimetry

[~, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

% Log tranform, zero mean, and normalise by SD
lsri = log(lsri);
lsri = lsri - mean(lsri(:,:),2);
lsri = lsri./std(lsri(:,:),[],2);

%% Perform CC algos

Lum = zeros(size(T_SRF,2),size(T_SPD,2));
for i=1:size(lsri,3)
    Lum(:,i) = T_lum'*(T_SRF.*T_SPD(:,i));
end

[~,KMM,~,~,~] = performCCandMarkForRange(lsri,Lum,1,100);

end