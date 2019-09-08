clear, clc, close all

plot_where = [500,200];
plot_size  = [800,400];

set(0,'defaultAxesFontName', 'Courier')

base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\melcomp_2_caller';

print_figures = 0;

%%
figure('Position',[plot_where plot_size]), hold on
melcomp_2('SPD','Granada_sub','plt','3D');
xlim([0.5 1])
cleanTicks
view(-170,40)

if print_figures
    save2pdf([base,'\viewpoint.pdf'])
end

%%

[~, T_SRF, T_SSF, ~, S_sh] = melcomp_loader('SRF','Vrhel_nat_2');

figure('Position',[plot_where plot_size]), hold on
%plot(SToWls(S_sh),T_SRF./T_SRF(12,:))
plot(SToWls(S_sh),T_SRF)
axis tight
xlim([min(xlim) 730])
%ylim([0 10])
%plot(SToWls(S_sh),T_SSF(:,[3,5])*max(ylim),'k-.')
cleanTicks
xticks('auto')
xlabel('Wavelength (nm)')
ylabel({'Reflectance or', 'relative spectral sensitivity'})

if print_figures
    save2pdf([base,'\plateau.pdf'])
end
