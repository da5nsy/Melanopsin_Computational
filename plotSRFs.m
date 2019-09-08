clear, clc, close all


%% Load Data

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_1',...
    'SSF','SS10',...
    'lum','CIE_10',...
    'mel_offset',0);

%%

figure, hold on
%plot(SToWls(S_sh),T_SRF./T_SRF(12,:))
plot(SToWls(S_sh),T_SRF)
axis tight
xlim([min(xlim) 730])
%ylim([0 10])
%plot(SToWls(S_sh),T_SSF(:,[3,5])*max(ylim),'k-.')
cleanTicks
xticks('auto')
xlabel('Wavelength (nm)')
ylabel({'SRF'})

%save2pdf('plateau.pdf')
