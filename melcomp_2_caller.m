clear, clc, close all

plot_where = [500,200];
plot_size  = [800,400];

set(0,'defaultAxesFontName', 'Courier')

base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\melcomp_2_caller';

print_figures = 1;

%%
figure('Position',[plot_where plot_size]), hold on
melcomp_2('SPD','Granada_sub','plt','3D');
xlim([0 1])
cleanTicks
view(-176,40)

if print_figures
    save2pdf([base,'\viewpoint.pdf'])
end
