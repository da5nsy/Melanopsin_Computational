%% 
clear, clc, close all

load T_melanopsin.mat

%% extracted from melcomp_6

figure(3)
cla
clear ex pc %only needed during debugging when rerunning script

pca_range = -70:5:130;

for i=1:length(pca_range)
    pc(i) = melcomp_6_looper_branch(pca_range(i),1,0);
    disp(pca_range(i))
end

[~, mel_peak_loc] = max(T_melanopsin);
S_melanopsin_f = SToWls(S_melanopsin);
mel_peak = S_melanopsin_f(mel_peak_loc);

for i=1:length(pca_range)
    ex(i) = pc(i).explained(3);
end

%figure, hold on
plot(pca_range+mel_peak,ex/max(ex),'k','DisplayName','PC3 score')
%plot
 
xlim('auto')
ylim([0 1])
yticks([min(ylim),max(ylim)])

legend('off')

xlabel('Wavelength shift (nm)')
ylabel('PC3 score')

%%

out = [pca_range;ex];

save('opt_melcomp6result_oldLMSToMacBoyn_oldrefs_sp_norm1','out')

