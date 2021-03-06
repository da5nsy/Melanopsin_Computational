% A fresh attempt at the optimality question

clear, clc, close all

mel_offset_range = [-100:5:200];
n=1;

for offset = mel_offset_range

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_1',...
    'SSF','SS10',...
    'mel_offset',offset,...
    'lum','CIE_10');

[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

pltc_alt = repmat(jet(size(T_SRF,2))',1,1,size(T_SPD,2));
rng(7); pltc_alt=pltc_alt(:,randperm(size(T_SRF,2)),:);

lsri = log(lsri);
lsri_c = lsri(:,:); 
for i=1:size(lsri,1)
    lsri_c(i,:) = (lsri_c(i,:) - mean(lsri_c(i,:)))./std(lsri_c(i,:));
end
lsri_c = reshape(lsri_c,size(lsri));


%[sf_l,sf_s] = melcomp_6_calcsf(lsri_c,0:0.01:1,-2:0.01:-0.5,1,pltc_alt); 
%[sf_l,sf_s] = melcomp_6_calcsf(lsri_c,0:0.01:1,-2:0.01:-0.5,0,pltc_alt);

[sf_l(n),sf_s(n)] = melcomp_6_calcsf(lsri_c,-5:0.01:5,-10:0.01:10,0,pltc_alt);

lsri_mel(:,:,:,n) = [lsri_c(1,:,:)+sf_l(n)*lsri_c(4,:,:);lsri_c(2,:,:)+sf_s(n)*lsri_c(4,:,:)];

K(n) = KMeansMark(lsri_mel(:,:,:,n));

n=n+1;

disp(offset)

end

%%
figure, hold on
plot(mel_offset_range +488,K)
[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_1',...
    'SSF','SS10',...
    'mel_offset',0,...
    'lum','CIE_10');
plot(SToWls(S_sh),T_SSF)

%%
figure, hold on

for i = 1:length(mel_offset_range)

lsri_mel_1 = lsri_mel(:,:,:,i);

cla
scatter(lsri_mel_1(1,:),lsri_mel_1(2,:),...
    [],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
axis equal
xlim([-10 10])
ylim([-10 10])
cleanTicks
title(mel_offset_range(i))

drawnow
pause(0.2)

end
